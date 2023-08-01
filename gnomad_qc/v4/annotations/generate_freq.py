"""
Script to generate the frequency data annotations across v4 exomes.

This script first splits the v4 VDS into multiples VDSs based which are then densified
and annotated with frequency data and histograms. The VDSs are then merged back together
in a hail Table. Next the script corrects for the high AB heterozygous GATK artifact in
existing annotations when given the AF threshold using a high AB het array. The script
then computes the inbreeding coefficient using the raw call stats. Finally, it computes
the filtering allele frequency and grpmax with the AB-adjusted frequencies.
"""
import argparse
import logging
from copy import deepcopy
from typing import List, Optional

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_downsamplings,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    merge_freq_arrays,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.filtering import split_vds_by_strata
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SORT_ORDER
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_freq
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency")
logger.setLevel(logging.INFO)

FREQ_ROW_FIELDS = [
    "freq",
    "high_ab_hets_by_group_membership",
    "qual_hists",
    "raw_qual_hists",
    "age_hist_het",
    "age_hist_hom",
]
"""
List of final top level row and global annotations created from dense data that we
want on the frequency HT before deciding on the AF cutoff.
"""

FREQ_GLOBAL_FIELDS = ["freq_meta", "age_distribution"]
"""
List of final global annotations created from dense data that we want on the frequency
HT before deciding on the AF cutoff.
"""


def get_freq_resources(
    overwrite: bool = False, test: Optional[bool] = False, chrom: Optional[str] = None
) -> PipelineResourceCollection:
    """
    Get frequency resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :param chrom: Chromosome used in freq calculations.
    :return: Frequency resources.
    """
    freq_pipeline = PipelineResourceCollection(
        pipeline_name="frequency",
        overwrite=overwrite,
    )
    run_freq_and_dense_annotations = PipelineStepResourceCollection(
        "--run-freq-and-dense-annotations",
        output_resources={
            "freq_and_dense_annotations": get_freq(
                test=test, hom_alt_adjusted=False, chrom=chrom
            ),
        },
    )
    correct_for_high_ab_hets = PipelineStepResourceCollection(
        "--correct-for-high-ab-hets",
        pipeline_input_steps=[run_freq_and_dense_annotations],
        output_resources={
            "freq_ht": get_freq(test=test, hom_alt_adjusted=True, chrom=chrom),
        },
    )
    freq_pipeline.add_steps(
        {
            "run_freq_and_dense_annotations": run_freq_and_dense_annotations,
            "correct_for_high_ab_hets": correct_for_high_ab_hets,
        }
    )
    return freq_pipeline


def get_vds_for_freq(
    use_test_dataset: hl.bool = False,
    test_gene: hl.bool = False,
    test_n_partitions: Optional[hl.int] = None,
    chrom: Optional[hl.int] = None,
) -> hl.vds.VariantDataset:
    """
    Prepare VDS for frequency calculation by filtering to release samples and only adding necessary annotations.

    :param use_test_dataset: Whether to use test dataset.
    :param test_gene: Whether to filter to DRD2 for testing purposes.
    :param test_n_partitions: Number of partitions to use for testing.
    :param chrom: Chromosome to filter to.
    :return: Hail VDS with only necessary annotations.
    """
    logger.info(
        "Reading the %s gnomAD v4 VDS...", "test" if use_test_dataset else "full"
    )
    if test_n_partitions:
        test_partitions = range(test_n_partitions)
    else:
        test_partitions = None

    vds = get_gnomad_v4_vds(
        test=use_test_dataset,
        release_only=True,
        filter_partitions=test_partitions,
        chrom=chrom,
        annotate_meta=True,
    )

    if test_gene:
        logger.info("Filtering to DRD2 in VDS for testing purposes...")
        test_interval = [
            hl.parse_locus_interval(
                "chr11:113409605-113475691", reference_genome="GRCh38"
            )
        ]
        vds = hl.vds.filter_intervals(vds, test_interval, split_reference_blocks=True)

    logger.info("Annotating VDS with only necessary sample metadata...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.select_cols(
            pop=vds.variant_data.meta.population_inference.pop,
            sex_karyotype=vds.variant_data.meta.sex_imputation.sex_karyotype,
            fixed_homalt_model=vds.variant_data.meta.project_meta.fixed_homalt_model,
            gatk_version=vds.variant_data.meta.project_meta.gatk_version,
            age=vds.variant_data.meta.project_meta.age,
            sample_age_bin=get_sample_age_bin(vds.variant_data.meta.project_meta.age),
            ukb_sample=hl.if_else(
                vds.variant_data.meta.project_meta.ukb_sample, "ukb", "non_ukb"
            ),
        ),
    )

    # Downsamplings are done outside the above function as we need to annotate
    # globals and rows.
    logger.info("Annotating downsampling groups...")
    vds.variant_data = annotate_downsamplings(
        vds.variant_data, DOWNSAMPLINGS["v4"], pop_expr=vds.variant_data.pop
    )

    logger.info("Annotating non_ref hets pre-split...")
    vds.variant_data = vds.variant_data.annotate_entries(
        _het_non_ref=vds.variant_data.LGT.is_het_non_ref()
    )

    logger.info("Selecting only required fields to reduce memory usage...")
    vds.variant_data = vds.variant_data.select_entries(
        "LA", "LAD", "DP", "GQ", "LGT", "_het_non_ref"
    )

    logger.info("Spltting mutliallelics in VDS...")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info(
        "Computing adj and _het_AD as part of reducing fields to reduce memory"
        " usage during dense dependent steps..."
    )
    vds = annotate_adj_and_select_fields(vds)

    return vds


def annotate_adj_and_select_fields(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Annotate adj, _het_ad, and select fields to reduce memory usage.

    :param vds: Hail VDS to annotate adj onto variant data.
    :return: Hail VDS with adj annotation.
    """
    rmt = vds.reference_data
    vmt = vds.variant_data

    rmt = rmt.annotate_entries(
        adj=(rmt.DP >= 10) & (rmt.GQ >= 20)
    )  # TODO: confirm this doesn't mess up haploids
    vmt = vmt.select_entries(
        "_het_non_ref",
        "DP",
        "GQ",
        "GT",
        adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD),
        _het_ad=vmt.AD[1],
    )
    return hl.vds.VariantDataset(rmt, vmt)


def correct_call_stats(ht: hl.Table, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct frequencies at sites with an AF greater than the af_threshold.

    :param ht: Hail Table containing freq and high_ab_het annotations.
    :param af_threshold: AF threshold at which to correct frequency. Default is 0.01.
    :return: Hail Table with adjusted frequencies.
    """
    ht = ht.annotate(
        ab_adjusted_freq=hl.if_else(
            ht.freq[0].AF > af_threshold,
            hl.map(
                lambda f, g: hl.struct(
                    AC=hl.int32(f.AC + g),
                    AN=f.AN,
                    homozygote_count=f.homozygote_count + g,
                    AF=hl.if_else(f.AN > 0, (f.AC + g) / f.AN, hl.missing(hl.tfloat64)),
                ),
                ht.freq,
                ht.high_ab_hets_by_group_membership,
            ),
            ht.freq,
        )
    )

    return ht


def create_high_ab_age_hists_expr(ht: hl.Table, age_group_key="sample_age_bin"):
    """
    Create histograms of high ab counts using age bins to account for high AB hets becoming hom alts.

    :param ht: Hail Table containing age hists, AB annotation.
    :param age_group_key: Age group key to use for age histogram.
    :return: Hail struct containing age histogram of high ab counts.
    """
    non_range_entries = hl.set(["n_larger", "n_smaller"])
    age_bins_indices = hl.sorted(
        hl.enumerate(ht["freq_meta"], index_first=False)
        .filter(lambda x: x[0].contains(age_group_key))
        .map(lambda x: (x[0][age_group_key], x[1]))
    )

    age_bins_indices_dict = hl.dict(age_bins_indices)
    age_bin_indices_no_edges = age_bins_indices.filter(
        lambda x: ~non_range_entries.contains(x[0])
    )
    return hl.struct(
        bin_freq=hl.starmap(
            lambda x, y: ht.high_ab_hets_by_group_membership[y],
            age_bin_indices_no_edges,
        ),
        n_smaller=ht.high_ab_hets_by_group_membership[
            age_bins_indices_dict["n_smaller"]
        ],
        n_larger=ht.high_ab_hets_by_group_membership[age_bins_indices_dict["n_larger"]],
    )


def get_sample_age_bin(
    sample_age: hl.expr.Int32Expression,
    lower_bound: int = 30,
    upper_bound: int = 80,
    number_of_bins: int = 10,
) -> hl.expr.StringExpression:
    """
    Get the age bin for a sample.

    :param sample_age: Sample age.
    :return: Sample age bin.
    """
    bin_size = (upper_bound - lower_bound) / number_of_bins
    lower_bin = hl.int(
        hl.floor((sample_age - lower_bound) / bin_size) * bin_size + lower_bound
    )
    upper_bin = hl.int(lower_bin + bin_size)

    bin_label = hl.if_else(
        sample_age < lower_bound,
        "n_smaller",
        hl.if_else(
            sample_age >= upper_bound,
            "n_larger",
            hl.str(lower_bin) + "-" + hl.str(upper_bin),
        ),
    )

    return hl.or_missing(hl.is_defined(sample_age), bin_label)


def compute_age_hist(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute age histograms for each variant.

    :param mt: Input MT with age annotation.
    :return: MatrixTable with age histogram annotations.
    """
    mt = mt.annotate_rows(**age_hists_expr(mt.adj, mt.GT, mt.age))

    # Compute callset-wide age histogram global
    mt = mt.annotate_globals(
        age_distribution=mt.aggregate_cols(hl.agg.hist(mt.age, 30, 80, 10))
    )
    return mt


def correct_age_hists(ht: hl.Table) -> hl.Table:
    """
    Correct age histograms.

    Correct by subtracting age_high_ab_hists from age_hist_het and adding
    age_high_ab_hists to age_hist_hom to account for high AB hets becoming hom alts.

    :param ht: Hail Table containing age hists and hist of AB counts by age annotation.
    :return: Hail Table
    """
    ht = ht.annotate(age_high_ab_hist=create_high_ab_age_hists_expr(ht))

    return ht.annotate(
        age_hist_het=hl.struct(
            bin_edges=ht.age_hist_het.bin_edges,
            bin_freq=hl.map(
                lambda x, y: x - y,
                ht.age_hist_het.bin_freq,
                ht.age_high_ab_hist.bin_freq,
            ),
            n_smaller=ht.age_hist_het.n_smaller - ht.age_high_ab_hist.n_smaller,
            n_larger=ht.age_hist_het.n_larger - ht.age_high_ab_hist.n_larger,
        ),
        age_hist_hom=hl.struct(
            bin_edges=ht.age_hist_hom.bin_edges,
            bin_freq=hl.map(
                lambda x, y: x + y,
                ht.age_hist_hom.bin_freq,
                ht.age_high_ab_hist.bin_freq,
            ),
            n_smaller=ht.age_hist_hom.n_smaller + ht.age_high_ab_hist.n_smaller,
            n_larger=ht.age_hist_hom.n_larger + ht.age_high_ab_hist.n_larger,
        ),
    )


def compute_qual_hists(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate quality metrics histograms.

    :param mt: Input MT.
    :return: MatrixTable with qual histogram annotations.
    """
    mt = mt.annotate_rows(
        qual_hists=qual_hist_expr(
            gt_expr=mt.GT,
            gq_expr=mt.GQ,
            dp_expr=mt.DP,
            adj_expr=mt.adj,
            ab_expr=mt._het_ad / mt.DP,
        )
    )
    mt = mt.annotate_rows(
        qual_hists=hl.Struct(
            **{
                i.replace("_adj", ""): mt.qual_hists[i]
                for i in mt.qual_hists
                if "_adj" in i
            }
        ),
        raw_qual_hists=hl.Struct(
            **{i: mt.qual_hists[i] for i in mt.qual_hists if "_adj" not in i}
        ),
    )
    return mt


def correct_qual_hists(ht: hl.Table) -> hl.Table:  # add ab_threshold as arg
    """
    Correct quality metrics histograms.

    Correct by accessing the qual_hist and raw_qual_hist structs and removing
    all counts from the ab_hist_alt array where bin_edges exceed 0.9 AB.

    :param ht: Hail Table containing qual hists, AB annotation.
    :return: Hail Table
    """

    def _correct_ab_hist_alt(ab_hist_alt):
        return hl.struct(
            bin_edges=ab_hist_alt.bin_edges,
            bin_freq=hl.map(
                lambda edge, freq: hl.if_else(edge >= 0.9, 0, freq),
                ab_hist_alt.bin_edges[:-1],
                ab_hist_alt.bin_freq,
            ),
            n_smaller=ab_hist_alt.n_smaller,
            n_larger=0,
        )

    qual_hists = ["qual_hists", "raw_qual_hists"]
    ht = ht.annotate(
        **{
            x: ht[x].annotate(ab_hist_alt=_correct_ab_hist_alt(ht[x].ab_hist_alt))
            for x in qual_hists
        }
    )
    return ht


def generate_faf_grpmax(ht: hl.Table) -> hl.Table:
    """
    Compute filtering allele frequencies and grpmax with the AB-adjusted frequencies.

    :param ht: Hail Table containing freq, ab_adjusted_freq, high_ab_het annotations.
    :return: Hail Table with faf & grpmax annotations.
    """
    faf, faf_meta = faf_expr(
        ht.ab_adjusted_freq, ht.freq_meta, ht.locus, POPS_TO_REMOVE_FOR_POPMAX
    )
    ht = ht.annotate(
        faf=faf,
        grpmax=pop_max_expr(
            ht.ab_adjusted_freq, ht.freq_meta, POPS_TO_REMOVE_FOR_POPMAX
        ),
    )
    ht = ht.annotate_globals(
        faf_meta=faf_meta,
        faf_index_dict=make_faf_index_dict(faf_meta, label_delimiter="-"),
    )
    ht = ht.annotate(
        grpmax=ht.grpmax.annotate(
            faf95=ht.faf[
                ht.faf_meta.index(lambda x: x.values() == ["adj", ht.grpmax.pop])
            ].faf95
        )
    )
    return ht


def generate_freq_and_hists_ht(
    vds: hl.vds.VariantDataset,
    ab_cutoff: float = 0.9,
    idx: Optional[int] = None,
) -> hl.Table:
    """
    Generate frequency and histogram annotations.

    Assumes all necessary annotations are present:
        - adj
        - _het_ad
        - _het_non_ref
        - GT
        - GQ
        - s
        - pop
        - sex_karyotype
        - fixed_homalt_model
        - gatk_version
        - age
        - sample_age_bin
        - ukb_sample
        - downsampling
        - downsamplings

    :param vds: Input VDS.
    :param ab_cutoff: Allele balance cutoff to use for high AB het annotation.
    :param idx: Optional index to append to temp file name.
    :return: Hail Table with frequency and histogram annotations.
    """
    final_rows_anns = {}
    final_globals_anns = {}

    logger.info("Densifying VDS # %s...", idx)
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Computing sex adjusted genotypes...")
    mt = mt.transmute_entries(
        GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype),
    )

    logger.info("Annotating frequencies and counting high AB het calls...")
    additional_strata_expr = [
        {"gatk_version": mt.gatk_version},
        {"gatk_version": mt.gatk_version, "pop": mt.pop},
        # TODO: We need this to get high AB hets for age histogram correction, dont
        #  care about the frequency of it, should we drop it?
        {"sample_age_bin": mt.sample_age_bin},
        # TODO: confirm this is how we want to name this.
        {"ukb_sample": mt.ukb_sample},
    ]

    def _needs_high_ab_het_fix(entry, col):
        return hl.int(
            entry.GT.is_het_ref()
            & (entry._het_ad / entry.DP > ab_cutoff)
            & entry.adj
            & ~col.fixed_homalt_model
            & ~entry._het_non_ref
        )  # Skip adjusting genotypes if sample originally had a het nonref genotype

    freq_ht = annotate_freq(
        mt,
        sex_expr=mt.sex_karyotype,
        pop_expr=mt.pop,
        downsamplings=hl.eval(mt.downsamplings),
        downsampling_expr=mt.downsampling,
        ds_pop_counts=hl.eval(mt.ds_pop_counts),
        additional_strata_expr=additional_strata_expr,
        entry_agg_funcs={
            "high_ab_hets_by_group_membership": (_needs_high_ab_het_fix, hl.agg.sum)
        },
        annotate_mt=False,
    )

    logger.info("Making freq index dict...")
    # Add additional strata to the sort order, keeping group, i.e. adj, at the end.
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version", "ukb_sample", "sample_age_bin"]

    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=freq_ht.freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
            # TODO: Check if we actually want to see age_bin, I dont think we do.
        )
    )
    logger.info("Setting Y metrics to NA for XX groups...")
    freq_ht = freq_ht.annotate(freq=set_female_y_metrics_to_na_expr(freq_ht))

    logger.info("Computing quality metrics histograms...")
    mt = compute_qual_hists(mt)

    logger.info("Computing age histograms for each variant...")
    mt = compute_age_hist(mt)  # global age distribution is a global

    hists = mt.rows()[freq_ht.key]
    final_rows_anns.update(
        {
            "qual_hists": hists.qual_hists,
            "raw_qual_hists": hists.raw_qual_hists,
            "age_hist_het": hists.age_hist_het,
            "age_hist_hom": hists.age_hist_hom,
        }
    )
    final_globals_anns.update({"age_distribution": mt.index_globals().age_distribution})

    freq_ht = freq_ht.annotate(**final_rows_anns)
    freq_ht = freq_ht.annotate_globals(**final_globals_anns)

    # TODO: Remove after testing.
    freq_ht.describe()
    freq_ht = freq_ht.checkpoint(
        new_temp_file(f"freq_ht_{idx}", extension="ht"),
        # f"gs://gnomad-mwilson/v4/frequencies/test/freq_ht_{idx}.ht",
        overwrite=args.overwrite,
        _read_if_exists=True,
    )

    return freq_ht


def merge_histograms(ht: hl.Table, indices: List[int]) -> hl.Table:
    """
    Merge histogram annotations.

    This function merges all split histogram annotations by
    summing the arrays in an element-wise fashion across the like histograms. Histograms
    that should be merged are named similarly with the index as a suffix. It keeps one
    bin_edge annotation but merges the bin_freq, n_smaller, and n_larger annotations by
    summing them.
    :param ht: Hail Table with histogram annotations.
    :param indices: List of indices to merge.
    :return: Hail Table with merged histogram annotations.
    """
    age_hists = ["age_hist_het", "age_hist_hom"]
    qual_hists = [
        "gq_hist_all",
        "dp_hist_all",
        "gq_hist_alt",
        "dp_hist_alt",
        "ab_hist_alt",
    ]
    hist_structs = {"qual_hists": qual_hists, "raw_qual_hists": qual_hists}

    def _hist_merge(arrays: List[hl.expr.StructExpression]):
        """
        Merge histograms.

        :param arrays: List of histogram structs to merge.
        :return: Merged histogram struct.
        """
        return hl.fold(
            lambda i, j: hl.struct(
                **{
                    "bin_edges": (
                        i.bin_edges
                    ),  # Bin edges are the same for all histograms
                    "bin_freq": hl.zip(i.bin_freq, j.bin_freq).map(
                        lambda x: x[0] + x[1]
                    ),
                    "n_smaller": i.n_smaller + j.n_smaller,
                    "n_larger": i.n_larger + j.n_larger,
                }
            ),
            arrays[0].select("bin_edges", "bin_freq", "n_smaller", "n_larger"),
            arrays[1:],
        )

    ht = ht.annotate(
        **{
            age_hist: _hist_merge([ht[f"{age_hist}_{idx}"] for idx in indices])
            for age_hist in age_hists
        },
        **{
            hist_struct: hl.struct(
                **{
                    hist: _hist_merge(
                        [ht[f"{hist_struct}_{idx}"][hist] for idx in indices]
                    )
                    for hist in hists
                }
            )
            for hist_struct, hists in hist_structs.items()
        },
    )
    ht = ht.annotate_globals(
        age_distribution=_hist_merge(
            [ht.index_globals()[f"age_distribution_{i}"] for i in indices]
        )
    )

    return ht


def combine_freq_hts(
    freq_hts: hl.DictExpression,
    row_annotations: List[str],
    globals_annotations: List[str],
) -> hl.Table:
    """
    Combine frequency HTs into a single HT.

    :param freq_hts: Dictionary of frequency HTs.
    :param row_annotations: List of annotations to put onto one hail Table.
    :param globals_annotations: List of global annotations to put onto one hail Table.
    :return: HT with all freq_hts annotations.
    """
    # Create new HT with just variants and downsamplings global annotation to join onto
    freq_ht = freq_hts[1].select().select_globals("downsamplings")

    # Annotate all hts' annotations with a idx suffix to the new HT. We don't remove
    # variants when splitting the VDS so each table has all rows.
    logger.info("Annotating frequency HT with split HT dense dependent annotations")
    freq_ht = freq_ht.annotate(
        **{
            f"{ann}_{idx}": freq_hts[idx][freq_ht.key][ann]
            for idx in freq_hts.keys()
            for ann in row_annotations
        }
    )
    freq_ht = freq_ht.annotate_globals(
        **{
            f"{ann}_{idx}": freq_hts[idx].index_globals()[ann]
            for idx in freq_hts.keys()
            for ann in globals_annotations
        }
    )

    # Combine freq arrays and high ab het counts by group arrays into single annotations
    logger.info(
        "Merging frequency arrays, metadata, and high ab het counts by group array..."
    )
    comb_freq, comb_freq_meta, comb_high_ab_hets = merge_freq_arrays(
        farrays=[freq_ht[f"freq_{idx}"] for idx in freq_hts.keys()],
        fmeta=[freq_ht[f"freq_meta_{idx}"] for idx in freq_hts.keys()],
        count_arrays=[
            freq_ht[f"high_ab_hets_by_group_membership_{idx}"]
            for idx in freq_hts.keys()
        ],
    )
    # TODO: We want to drop all _idx annotations here but keeping them around for now for testing
    # This could also live in the merge function, passed as a boolean parameter
    freq_ht = freq_ht.annotate(
        freq=comb_freq,
        high_ab_hets_by_group_membership=comb_high_ab_hets,
    )

    # Merge all histograms into single annotations
    logger.info("Merging all histograms...")
    freq_ht = merge_histograms(freq_ht, freq_hts.keys())

    logger.info("Final frequency HT schema...")
    freq_ht.describe()
    logger.info("Making freq index dict...")
    # Add our additional strata to the sort order, keeping group, i.e. adj, at the end
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version", "ukb_sample", "sample_age_bin"]
    freq_ht = freq_ht.annotate_globals(freq_meta=hl.eval(comb_freq_meta))
    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=freq_ht.freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
            # TODO: Check if we actually want to see age_bin, I dont thin we do
        ),
    )

    return freq_ht


def main(args):  # noqa: D103
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = use_test_dataset or test_n_partitions or test_gene
    chrom = args.chrom
    ab_cutoff = args.ab_cutoff
    af_threshold = args.af_threshold
    correct_for_high_ab_hets = args.correct_for_high_ab_hets

    hl.init(
        log="/generate_frequency_data.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    resources = get_freq_resources(args.overwrite, test, chrom)

    if args.run_freq_and_dense_annotations:
        logger.info("Running dense dependent steps...")
        res = resources.run_freq_and_dense_annotations
        res.check_resource_existence()

        logger.info(
            "Getting multi-allelic split VDS with adj and _het_AD entry annotations..."
        )
        vds = get_vds_for_freq(use_test_dataset, test_gene, test_n_partitions, chrom)

        if args.split_vds_by_annotation:
            logger.info(
                "Splitting VDS by ukb_sample annotation to reduce data size for"
                " densification..."
            )
            vds_list = split_vds_by_strata(vds, strata_expr=vds.variant_data.ukb_sample)
            freq_hts = {}

            # Make a dict of freq hts with idx as key which will be used when joining
            # annotations. This is more straight forward to me than using joins since
            # hail keeps the existing annotation on the left table and the right table
            # gets a suffix added to the same annotation. This way, every table's
            # annotation will end up with a suffix in the merge.
            for idx, vds in enumerate(vds_list, start=1):
                freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff, idx=idx)
                freq_hts[idx] = freq_ht

            freq_ht = combine_freq_hts(freq_hts, FREQ_ROW_FIELDS, FREQ_GLOBAL_FIELDS)
        else:
            freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff)

        freq_ht.write(res.freq_and_dense_annotations.path, overwrite=args.overwrite)

    if correct_for_high_ab_hets:
        logger.info(
            "Adjusting annotations impacted by high AB het -> hom alt adjustment..."
        )
        res = resources.correct_for_high_ab_hets
        res.check_resource_existence()
        ht = res.freq_and_dense_annotations.ht()

        logger.info("Correcting call stats...")
        ht = correct_call_stats(ht, af_threshold)

        logger.info("Correcting qual AB histograms...")
        ht = correct_qual_hists(ht)

        logger.info("Correcting age histograms...")
        ht = correct_age_hists(ht)

        logger.info("computing FAF & grpmax...")
        ht = generate_faf_grpmax(ht)

        logger.info("Calculating InbreedingCoeff...")
        ht = ht.annotate(
            InbreedingCoeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1])
        )

        # TODO: Leaving in know while we test but need to drop fields we do not want
        # -- 'age_high_ab_his', all annotations from the split VDSs, only keep combinged,
        # rename ab_adjusted_freq to just freq and decide if we want to store uncorrect?
        # Probably,just rename it, also remove age bins and gatk versions from freq fields?
        # Also, change captialization of hists depending on decision from DP slack
        # thread
        logger.info("Writing frequency table...")
        ht.describe()
        ht.write(res.freq_ht.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--use-test-dataset",
        help="Runs a test on the gnomad test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--test-gene",
        help="Runs a test on the DRD2 gene in the gnomad test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the VDS as input for testing purposes. Defaults"
            "to 2 if passed without a value."
        ),
        nargs="?",
        const=2,
        type=int,
    )
    parser.add_argument(
        "--chrom",
        help="If passed, script will only run on passed chromosome.",
        type=str,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--run-freq-and-dense-annotations",
        help=(
            "Calculate frequencies, histograms, and high AB sites per sample grouping."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--split-vds-by-annotation",
        help=(
            "Split VDS by annotation to reduce data size for densification."
            " Defaults to False."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--correct-for-high-ab-hets",
        help=(
            "Correct each frequency entry to account for homozygous alternate depletion"
            " present in GATK versions released prior to 4.1.4.1 and run chosen"
            " downstream annotations."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--ab-cutoff",
        help=(
            "Allele balance threshold to use when adjusting heterozygous calls to "
            "homozygous alternate calls at sites for samples that used GATK versions"
            " released prior to 4.1.4.1."
        ),
        type=float,
        default=0.9,
    )
    parser.add_argument(
        "--af-threshold",
        help=(
            "Threshold at which to adjust site group frequencies at sites for"
            " homozygous alternate depletion present in GATK versions released prior to"
            " 4.1.4.1."
        ),
        type=float,
        default=0.01,
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
