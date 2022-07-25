import argparse
import functools
import logging
import operator
from typing import Callable, List

import hail as hl

from gnomad.sample_qc.pipeline import annotate_sex
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
    ukb_f_stat,
)
from gnomad_qc.v4.resources.sample_qc import (
    f_stat_sites,
    interval_coverage,
    platform,
    sex,
    sex_imputation_coverage,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


def determine_fstat_sites(
    vds: hl.vds.VariantDataset,
    approx_af_and_no_callrate: bool = False,
    min_af: float = 0.001,
    min_callrate: float = 0.99,
) -> hl.Table:
    """
    Write a Table with chromosome X SNPs that are bi-allelic, common, and high callrate by default.

    This Table is designed to be used as a variant filter in sex imputation for f-stat computation.

    .. warning::

        By default `approx_af_and_no_callrate` is False and the final Table will be filtered to high callrate (> value
        specified by `min_callrate`) variants. This requires a densify of chrX!!"

    .. note::

        If `approx_af_and_no_callrate` is True allele frequency is approximated with AC/(n_samples * 2) and no callrate
        filter is used.

    :param vds: Input VariantDataset
    :param approx_af_and_no_callrate: Whether to approximate allele frequency with AC/(n_samples * 2) and use no
        callrate cutoff to filter sites
    :param min_af: Alternate allele frequency cutoff used to filter sites
    :param min_callrate: Callrate cutoff used to filter sites
    :return: Table of chromosome X sites to be used for f-stat computation
    """

    vds = hl.vds.filter_chromosomes(vds, keep=["chrX"])
    vd = vds.variant_data
    vd = vd.filter_rows(
        (hl.len(vd.alleles) == 2) & hl.is_snp(vd.alleles[0], vd.alleles[1])
    )
    vd = vd.transmute_entries(GT=hl.experimental.lgt_to_gt(vd.LGT, vd.LA))
    vds = hl.vds.VariantDataset(vds.reference_data, vd)

    if approx_af_and_no_callrate:
        n_samples = vd.count_cols()
        logger.info("Number of samples: %d", n_samples)
        ht = vd.select_rows(
            AF=hl.agg.sum(vd.GT.n_alt_alleles()) / (n_samples * 2),
        ).rows()
        ht = ht.filter(ht.AF > min_af)
        ht = ht.annotate_globals(
            min_approx_af=min_af,
        )
    else:
        mt = hl.vds.to_dense_mt(vds)
        ht = hl.variant_qc(mt).rows()
        ht = ht.filter(
            (ht.variant_qc.call_rate > min_callrate) & (ht.variant_qc.AF[1] > min_af)
        )
        ht = ht.annotate(AF=ht.variant_qc.AF[1])
        ht = ht.annotate_globals(
            min_af=min_af,
            min_callrate=min_callrate,
        )

    return ht


def load_platform_ht(
    test: bool = False,
    calling_interval_name: str = "intersection",
    calling_interval_padding: int = 50,
) -> hl.Table:
    """
    Load platform assignment Table or test Table and return an error if requested Table does not exist.

    :param test: Whether a test platform assignment Table should be loaded if the final platform assignment Table doesn't exist
    :param calling_interval_name: Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or
        'intersection'. Only used if `test` is True and the final platform assignment Table is not already written
    :param calling_interval_padding: Number of base pair padding to use on the calling intervals. One of 0 or 50 bp.
        Only used if `test` is True and the final platform assignment Table is not already written
    :return: Platform Table
    """
    logger.info("Loading platform information...")
    if file_exists(platform.path):
        ht = platform.ht()
    elif test:
        test_platform_path = get_checkpoint_path(
            f"test_platform_assignment.{calling_interval_name}.pad{calling_interval_padding}"
        )
        if file_exists(test_platform_path):
            ht = hl.read_table(test_platform_path)
        else:
            raise FileNotFoundError(
                f"There is no final platform assignment Table written and a test platform assignment Table "
                f"does not exist for calling interval {calling_interval_name} and interval padding "
                f"{calling_interval_padding}. Please run platform_inference.py --assign-platforms "
                f"with the --test argument and needed --calling-interval-name/--calling-interval-padding "
                f"arguments."
            )
    else:
        raise FileNotFoundError(
            f"There is no final platform assignment Table written. Please run: "
            f"platform_inference.py --assign-platforms to compute the platform assignment Table."
        )

    return ht


def load_coverage_mt(
    test: bool = False,
    calling_interval_name: str = "intersection",
    calling_interval_padding: int = 50,
) -> hl.Table:
    """
    Load interval coverage MatrixTable or test MatrixTable and return an error if requested MatrixTable does not exist.

    :param test: Whether a test interval coverage MatrixTable should be loaded if the final Table doesn't exist
    :param calling_interval_name: Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or
        'intersection'. Only used if `test` is True and the final coverage MatrixTable is not already written
    :param calling_interval_padding: Number of base pair padding to use on the calling intervals. One of 0 or 50 bp.
        Only used if `test` is True and the final coverage MatrixTable is not already written
    :return: Interval coverage Table
    """
    logger.info("Loading interval coverage MatrixTable...")
    if file_exists(interval_coverage.path):
        mt = interval_coverage.mt()
    elif test:
        test_coverage_path = get_checkpoint_path(
            f"test_interval_coverage.{calling_interval_name}.pad{calling_interval_padding}",
            mt=True,
        )
        if file_exists(test_coverage_path):
            mt = hl.read_matrix_table(test_coverage_path)
        else:
            raise FileNotFoundError(
                f"There is no final coverage MatrixTable written and a test interval coverage MatrixTable does "
                f"not exist for calling interval {calling_interval_name} and interval padding "
                f"{calling_interval_padding}. Please run hard_filters.py --compute-coverage with the "
                f"--test argument and needed --calling-interval-name/--calling-interval-padding arguments."
            )
    else:
        raise FileNotFoundError(
            f"There is no final coverage MatrixTable written. Please run: "
            f"hard_filters.py --compute-coverage to compute the interval coverage MatrixTable."
        )

    return mt


def generate_sex_imputation_interval_qc_mt(
    mt: hl.MatrixTable,
    platform_ht: hl.Table,
    contigs: List[str] = ["chrX", "chrY", "chr20"],
    mean_dp_thresholds: List[int] = [5, 10, 15, 20, 25],
) -> hl.Table:
    """
    Create a Table of fraction of samples per interval and per platform with mean DP over specified thresholds.

    :param mt: Input interval coverage MatrixTable
    :param platform_ht: Input platform assignment Table
    :param contigs: Which contigs to filter to before computing the fraction of sample over DP thresholds
    :param mean_dp_thresholds: List of mean DP thresholds to use for computing the fraction of samples with mean
        interval DP >= the threshold
    :return: Table with annotations for the fraction of samples per interval and per platform over DP thresholds
    """
    logger.info(
        "Filtering interval coverage MatrixTable to the following contigs: %s...",
        ", ".join(contigs),
    )
    contigs = hl.literal(contigs)
    mt = mt.filter_rows(contigs.contains(mt.interval.start.contig))
    mt = mt.annotate_cols(platform=platform_ht[mt.col_key].qc_platform)
    mt = mt.annotate_rows(
        **{
            f"over_{dp}x": hl.agg.fraction(mt.mean_dp >= dp)
            for dp in mean_dp_thresholds
        }
    )

    logger.info("Adding per platform aggregation...")
    mt = mt.group_cols_by(mt.platform).aggregate(
        **{
            f"platform_over_{dp}x": hl.agg.fraction(mt.mean_dp >= dp)
            for dp in mean_dp_thresholds
        }
    )

    platform_ht = platform_ht.group_by(platform_ht.qc_platform).aggregate(
        n_samples=hl.agg.count()
    )
    mt = mt.annotate_cols(n_samples=platform_ht[mt.col_key].n_samples)
    mt = mt.select_globals(
        "calling_interval_name",
        "calling_interval_padding",
        mean_dp_thresholds=mean_dp_thresholds,
    )

    return mt


def compute_sex(
    vds: hl.vds.VariantDataset,
    coverage_mt: hl.MatrixTable,
    high_cov_intervals: bool = False,
    per_platform: bool = False,
    high_cov_by_platform_all: bool = False,
    platform_ht: hl.Table = None,
    min_platform_size: bool = 100,
    normalization_contig: str = "chr20",
    variant_depth_only_x_ploidy: bool = False,
    variant_depth_only_y_ploidy: bool = False,
    x_cov: int = None,
    y_cov: int = None,
    norm_cov: int = None,
    prop_samples_x: float = None,
    prop_samples_y: float = None,
    prop_samples_norm: float = None,
    freq_ht=None,
    min_af: float = 0.001,
    f_stat_cutoff: float = 0.5,
) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Within this function there are some critical parameters to consider (more details on each are described below):
        - Filtering to intervals with high coverage for sex chromosome ploidy estimation (`high_cov_intervals`, `x_cov`,
         `y_cov`, `norm_cov`, `prop_samples_x`, `prop_samples_y`, and `prop_samples_norm`).
        - Per platform computations: either to determine per platform sex karyotype cutoffs, or if `high_cov_intervals`
         is used, per platform sex chromosome ploidy estimation and determination of sex karyotype cutoffs is performed
         using per platform high coverage intervals (`x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`, `prop_samples_y`,
         `prop_samples_norm`, and `min_platform_size`).
        - Use per platform stats to determine which are high coverage across all platforms (depending on platform sample
         size). Uses parameters `high_cov_by_platform_all`, `x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
         `prop_samples_y`, `prop_samples_norm`.
        - We use `annotate_sex`, which uses Hail's VDS imputation of sex ploidy `hail.vds.impute_sex_chromosome_ploidy`,
          and uses `variant_depth_only_x_ploidy` and `variant_depth_only_y_ploidy`.
        - Use of a `freq_ht` to filter variants for the X-chromosome heterozygosity computation. This is the f_stat
         annotation applied by Hail's `impute_sex` module. (`freq_ht`, `min_af`).

    The following are the available options for chrX and chrY relative sex ploidy estimation and sex karyotype
    cutoff determination:
        - Ploidy estimation using all intervals found in `coverage_mt`.
        - Per platform ploidy estimation using all intervals found in `coverage_mt`. Sex karyotype cutoffs are
        determined by per platform ploidy distributions.
        - Ploidy estimation using only high coverage intervals across the entire dataset, defined as: chrX intervals
         having a proportion of all samples (`prop_samples_x`) with over a specified mean DP (`x_cov`), chrY intervals
         having a proportion of all samples (`prop_samples_y`) with over a specified mean DP (`y_cov`), and
         `normalization_contig` intervals having a proportion of all samples (`prop_samples_norm`) with over a
         specified mean DP (`norm_cov`). Sex karyotype cutoffs are determined by ploidy distributions of all samples.
        - Per platform ploidy estimation using only per platform high coverage intervals defined by: chrX intervals with
         the specific platform having a proportion of samples (`prop_samples_x`) with over a specified mean DP (`x_cov`),
         chrY intervals with the specific platform having a proportion of samples (`prop_samples_y`) with over a
         specified mean DP (`y_cov`), and `normalization_contig` intervals with the specific platform having a
         proportion of samples (`prop_samples_norm`) with over a specified mean DP (`norm_cov`). Sex karyotype cutoffs
         are determined by per platform ploidy distributions.
        - A combined ploidy estimation using only per platform high coverage intervals. Each interval is assessed as high
         coverage per platform on chrX, chrY, and `normalization_contig` as described above for "per platform ploidy
         estimation". Then this method uses only platforms that have # of samples > `min_platform_size` to determine
         intervals that have a high coverage across all platforms. Then sex ploidy estimation and sex karyotype cutoffs
         are determined using this intersection of high coverage intervals across platforms.

    For each of the options described above, there is also the possibility to use a ploidy estimation that uses only
    variants within the specified calling intervals:
        - This can be defined differently for chrX and chrY using `variant_depth_only_x_ploidy` and
        `variant_depth_only_y_ploidy`.
        - `variant_depth_only_x_ploidy` is the preferred method for chrX ploidy computation.
        - If not using only variants Hail's `hail.vds.impute_sex_chromosome_ploidy` method will only compute chromosome
         ploidy using reference block DP per calling interval. This method breaks up the reference blocks at the
         calling interval boundries, maintaining all reference block END information for the mean DP per interval
         computation. This is different from the method used on sparse MatrixTables in gnomad_methods
         `gnomad.utils.sparse_mt.impute_sex_ploidy`. In this the chromosome coverage will be computed using all
         reference blocks and variants found in the sparse MatrixTable, but it only uses specified calling intervals to
         determine the contig size and doesn't break up the reference blocks in the same way the Hail method does.

    :param vds: Input VDS for use in sex inference
    :param coverage_mt: Input interval coverage MatrixTable
    :param high_cov_intervals: Whether to filter to high coverage intervals for the sex ploidy and karyotype inference
    :param per_platform: Whether to run the sex ploidy and karyotype inference per platform
    :param high_cov_by_platform_all: Whether to filter to high coverage intervals for the sex ploidy and karyotype
        inference. Using only intervals that are considered high coverage across all platforms
    :param platform_ht: Input platform assignment Table. This is only needed if per_platform or high_cov_by_platform_all are True
    :param min_platform_size: Required size of a platform to be considered when using `high_cov_by_platform_all`. Only
        platforms that have # of samples > 'min_platform_size' are used to determine intervals that have a
        high coverage across all platforms
    :param normalization_contig: Which autosomal chromosome to use for normalizing the coverage of chromosomes X and Y
    :param variant_depth_only_x_ploidy: Whether to use depth of variant data within calling intervals instead of reference
        data for chrX ploidy estimation. Default will only use reference data
    :param variant_depth_only_y_ploidy: Whether to use depth of variant data within calling intervals instead of reference
        data for chrY ploidy estimation. Default will only use reference data
    :param x_cov: Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the
        interval_coverage MatrixTable
    :param y_cov: Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the
        interval_coverage MatrixTable
    :param norm_cov: Mean coverage level used to define high coverage intervals on the normalization autosome
        (`normalization_contig`). This field must be in the interval_coverage MT
    :param prop_samples_x: Proportion samples at specified coverage `x_cov` to determine high coverage intervals on chromosome X
    :param prop_samples_y: Proportion samples at specified coverage `y_cov` to determine high coverage intervals on chromosome Y
    :param prop_samples_norm: Proportion samples at specified coverage `norm_cov` to determine high coverage intervals
        on the normalization chromosome specified by `normalization_contig`
    :param freq_ht: Table to use for f-stat allele frequency cutoff. The input VDS is filtered to sites in this Table
        prior to running Hail's `impute_sex` module, and alternate allele frequency is used from this Table with a
        `min_af` cutoff
    :param min_af: Minimum alternate allele frequency to be used in f-stat calculations
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff
    :return: Table with inferred sex annotation
    """
    if (per_platform or high_cov_by_platform_all) and platform_ht is None:
        raise ValueError(
            "'platform_ht' must be defined if 'per_platform' or 'high_cov_by_platform_all' is True!"
        )

    if per_platform and high_cov_by_platform_all:
        raise ValueError(
            "Only one of 'per_platform' or 'high_cov_by_platform_all' can be True!"
        )

    if high_cov_intervals and high_cov_by_platform_all:
        logger.warning(
            "Both 'high_cov_intervals' and 'high_cov_by_platform_all' are True, high coverage intervals will be defined"
            " as intervals that match high coverage criteria within all platforms."
        )

    def _get_high_coverage_intervals_ht(
        coverage_mt: hl.MatrixTable,
        prefix: str = "",
        agg_func: Callable[
            [hl.expr.BooleanExpression], hl.BooleanExpression
        ] = hl.agg.all,
    ) -> hl.Table:
        """
        Helper function to create a Table filtered to high coverage intervals.

        High coverage intervals are determined using `agg_func`, `x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
        `prop_samples_y`, and `prop_samples_norm`.

        :param coverage_mt: Input interval coverage MatrixTable
        :param prefix: Prefix of annotations in `coverage_mt` that contain the proportion of samples with mean DP over
            coverage cutoffs
        :param agg_func: Hail aggregation function to determine if an interval coverage meets the
            `cov_*` > `prop_samples_*` criteria
        :return: Table of high coverage intervals
        """
        prop_samples_dict = {
            "chrX": (x_cov, prop_samples_x),
            "chrY": (y_cov, prop_samples_y),
            normalization_contig: (norm_cov, prop_samples_norm),
        }
        filter_expr = functools.reduce(
            operator.or_,
            [
                (coverage_mt.interval.start.contig == contig)
                & agg_func(coverage_mt[f"{prefix}over_{cov}x"] > prop_samples)
                for contig, (cov, prop_samples) in prop_samples_dict.items()
            ],
        )

        return coverage_mt.filter_rows(filter_expr).rows()

    def _annotate_sex(vds, calling_intervals_ht):
        """
        Helper function to perform `annotate_sex` using unchanged parameters with changes to the VDS and calling intervals.

        :param vds: Input VDS to use for sex annotation
        :param calling_intervals_ht: Calling intervals to filter to for sex annotation
        :return: Table containing sex annotation for samples in the input VDS
        """
        return annotate_sex(
            vds,
            included_intervals=calling_intervals_ht,
            normalization_contig=normalization_contig,
            sites_ht=freq_ht,
            aaf_expr="AF",
            gt_expr="LGT",
            f_stat_cutoff=f_stat_cutoff,
            aaf_threshold=min_af,
            variants_only_x_ploidy=variant_depth_only_x_ploidy,
            variants_only_y_ploidy=variant_depth_only_y_ploidy,
        )

    if high_cov_intervals or high_cov_by_platform_all:
        logger.info(
            "Running sex ploidy and sex karyotype estimation using only high coverage intervals: %d percent of samples "
            "with greater than %dx coverage on chrX, %d percent of samples with greater than %dx coverage on chrY, and "
            "%d percent of samples with greater than %dx coverage on %s...",
            int(prop_samples_x * 100),
            x_cov,
            int(prop_samples_y * 100),
            y_cov,
            int(prop_samples_norm * 100),
            norm_cov,
            normalization_contig,
        )

    if per_platform:
        logger.info(
            "Running sex ploidy and sex karyotype estimation per platform %s...",
            "using per platform high coverage intervals" if high_cov_intervals else "",
        )
        platforms = platform_ht.aggregate(
            hl.agg.collect_as_set(platform_ht.qc_platform)
        )
        per_platform_sex_hts = []
        x_ploidy_cutoffs = {}
        y_ploidy_cutoffs = {}

        for platform in platforms:
            logger.info("Performing estimation for platform %s...", platform)
            if high_cov_intervals:
                calling_intervals_ht = _get_high_coverage_intervals_ht(
                    coverage_mt.filter_cols(coverage_mt.platform == platform),
                    prefix="platform_",
                )
            else:
                calling_intervals_ht = coverage_mt.rows()
            sex_ht = _annotate_sex(
                hl.vds.filter_samples(
                    vds, platform_ht.filter(platform_ht.qc_platform == platform)
                ),
                calling_intervals_ht,
            )
            sex_ht = sex_ht.annotate(platform=platform)
            per_platform_sex_hts.append(sex_ht)
            x_ploidy_cutoffs[platform] = sex_ht.index_globals().x_ploidy_cutoffs
            y_ploidy_cutoffs[platform] = sex_ht.index_globals().y_ploidy_cutoffs
        sex_ht = per_platform_sex_hts[0].union(*per_platform_sex_hts[1:])
        sex_ht = sex_ht.annotate_globals(
            x_ploidy_cutoffs=hl.struct(**x_ploidy_cutoffs),
            y_ploidy_cutoffs=hl.struct(**y_ploidy_cutoffs),
        )
    elif high_cov_by_platform_all:
        logger.info(
            "Running sex ploidy and sex karyotype estimation using high coverage intervals across all platforms. "
            "Limited to platforms with at least %s samples...",
            min_platform_size,
        )
        coverage_mt = coverage_mt.filter_cols(
            coverage_mt.n_samples >= min_platform_size
        )
        calling_intervals_ht = _get_high_coverage_intervals_ht(
            coverage_mt, prefix="platform_"
        )
        sex_ht = _annotate_sex(vds, calling_intervals_ht)
    elif high_cov_intervals:
        logger.info(
            "Running sex ploidy and sex karyotype estimation using high coverage intervals across the full sample set..."
        )
        calling_intervals_ht = _get_high_coverage_intervals_ht(
            coverage_mt, agg_func=lambda x: x
        )
        sex_ht = _annotate_sex(vds, calling_intervals_ht)
    else:
        logger.info("Running sex ploidy and sex karyotype estimation...")
        sex_ht = _annotate_sex(vds, coverage_mt.rows())

    sex_ht = sex_ht.annotate_globals(
        high_cov_intervals=high_cov_intervals,
        per_platform=per_platform,
        high_cov_by_platform_all=high_cov_by_platform_all,
        min_platform_size=min_platform_size,
        normalization_contig=normalization_contig,
        variant_depth_only_x_ploidy=variant_depth_only_x_ploidy,
        variant_depth_only_y_ploidy=variant_depth_only_y_ploidy,
        x_cov=x_cov,
        y_cov=y_cov,
        norm_cov=norm_cov,
        prop_samples_x=prop_samples_x,
        prop_samples_y=prop_samples_y,
        prop_samples_norm=prop_samples_norm,
        f_stat_min_af=min_af,
        f_stat_cutoff=f_stat_cutoff,
    )

    return sex_ht


def main(args):
    hl.init(log="/sex_inference.log", default_reference="GRCh38")

    test = args.test
    calling_interval_name = args.calling_interval_name
    calling_interval_padding = args.calling_interval_padding
    normalization_contig = args.normalization_contig

    try:
        if args.determine_fstat_sites:
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            ht = determine_fstat_sites(
                vds,
                approx_af_and_no_callrate=args.approx_af_and_no_callrate,
                min_af=args.min_af,
                min_callrate=args.min_callrate,
            )
            ht.naive_coalesce(args.fstats_n_partitions).write(
                get_checkpoint_path("test_f_stat_sites") if test else f_stat_sites.path,
                overwrite=args.overwrite,
            )

        if args.sex_imputation_interval_qc:
            coverage_mt = load_coverage_mt(
                test, calling_interval_name, calling_interval_padding
            )
            platform_ht = load_platform_ht(
                test, calling_interval_name, calling_interval_padding
            )
            coverage_mt = generate_sex_imputation_interval_qc_mt(
                coverage_mt,
                platform_ht,
                contigs=["chrX", "chrY", normalization_contig],
                mean_dp_thresholds=args.mean_dp_thresholds,
            )
            coverage_mt.naive_coalesce(args.interval_qc_n_partitions).write(
                get_checkpoint_path("test_sex_imputation_cov", mt=True)
                if test
                else sex_imputation_coverage.path,
                overwrite=args.overwrite,
            )

        if args.impute_sex:
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            platform_ht = load_platform_ht(
                test, calling_interval_name, calling_interval_padding
            )
            if args.f_stat_ukbb_var:
                # The UK Biobank f-stat table contains only variants that were high callrate (0.99) and common
                # (AF >0.001) within the UK Biobank 200K Regeneron exome dataset and it includes the UK Biobank 200K
                # allele frequency information that can be used in the hl.impute_sex f-stat computation allele frequency
                # cutoff (args.min-af)
                freq_ht = ukb_f_stat.ht()
            else:
                freq_ht = (
                    hl.read_table(get_checkpoint_path("test_f_stat_sites"))
                    if test
                    else f_stat_sites.ht()
                )

            sex_ht_path = (
                get_checkpoint_path(
                    f"sex_imputation{'.per_platform' if args.per_platform else ''}"
                    f"{'.high_cov_by_platform_all' if args.high_cov_by_platform_all else ''}"
                    f"{'.high_cov' if args.high_cov_intervals else ''}"
                    f"{'.ukbb_f_stat' if args.f_stat_ukbb_var else ''}"
                )
                if test
                else sex.path
            )

            # Added because without this impute_sex_chromosome_ploidy will still run even with overwrite=False
            if args.overwrite or not file_exists(sex_ht_path):
                coverage_mt = (
                    hl.read_matrix_table(
                        get_checkpoint_path("test_sex_imputation_cov", mt=True)
                    )
                    if test
                    else sex_imputation_coverage.mt()
                )

                ht = compute_sex(
                    vds,
                    coverage_mt,
                    args.high_cov_intervals,
                    args.per_platform,
                    args.high_cov_by_platform_all,
                    platform_ht,
                    args.min_platform_size,
                    normalization_contig,
                    args.variant_depth_only_x_ploidy,
                    args.variant_depth_only_y_ploidy,
                    args.x_cov,
                    args.y_cov,
                    args.norm_cov,
                    args.prop_samples_x,
                    args.prop_samples_y,
                    args.prop_samples_norm,
                    freq_ht,
                    args.min_af,
                    args.f_stat_cutoff,
                )
                ht = ht.annotate_globals(
                    f_stat_ukbb_var=args.f_stat_ukbb_var,
                )
                logger.info("Writing sex inference Table...")
                ht.write(sex_ht_path, overwrite=True)
            else:
                logger.warning("File exists and overwrite is not set!")
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("sex_inference"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Test the pipeline using the gnomAD v4 test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--determine-fstat-sites",
        help=(
            "Create Table of common (> value specified by '--min-af'), bi-allelic SNPs on chromosome X for f-stat "
            "calculations. Additionally filter to high callrate (> value specified by '--min-callrate') variants "
            "if '--approx-af-and-no-callrate' is not used. NOTE: This requires a densify of chrX!!"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--approx-af-and-no-callrate",
        help=(
            "Whether to approximate allele frequency with AC/(n_samples * 2) and use no callrate cutoff for "
            "determination of f-stat sites."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--fstat-n-partitions",
        help="Number of desired partitions for the f-stat sites output Table.",
        default=1000,
        type=int,
    )
    parser.add_argument(
        "--sex-imputation-interval-qc",
        help=(
            "Create a Table of the fraction of samples per interval and per platform with mean DP over thresholds "
            "specified by '--mean-dp-thresholds'."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--interval-qc-n-partitions",
        help="Number of desired partitions for the sex imputation interval QC output Table.",
        default=500,
        type=int,
    )
    parser.add_argument(
        "--normalization-contig",
        help="Which autosomal chromosome to use for normalizing the coverage of chromosomes X and Y.",
        type=str,
        default="chr20",
    )
    parser.add_argument(
        "--mean-dp-thresholds",
        help=(
            "List of mean DP cutoffs to determining the fraction of samples with mean coverage >= the cutoff for each "
            "interval."
        ),
        type=int,
        nargs="+",
        default=[5, 10, 15, 20, 25],
    )
    parser.add_argument(
        "--impute-sex",
        help="Run sex ploidy and sex karyotyping imputation.",
        action="store_true",
    )
    parser.add_argument(
        "--f-stat-ukbb-var",
        help=(
            "Whether to use UK Biobank high callrate (0.99) and common variants (UKBB allele frequency > value specified "
            "by '--min-af') for f-stat computation instead of the sites determined by '--determine-fstat-sites'. "
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-callrate", help="Minimum variant callrate.", default=0.99, type=float
    )
    parser.add_argument(
        "--min-af",
        help="Minimum variant allele frequency to retain variant in qc matrix table.",
        default=0.001,
        type=float,
    )
    parser.add_argument(
        "--f-stat-cutoff",
        help=(
            "Cutoff for f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY "
            "are above cutoff."
        ),
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--high-cov-intervals",
        help="Whether to filter to high coverage intervals for the sex ploidy and karyotype inference.",
        action="store_true",
    )
    parser.add_argument(
        "--per-platform",
        help="Whether to run the sex ploidy and karyotype inference per platform.",
        action="store_true",
    )
    parser.add_argument(
        "--high-cov-by-platform-all",
        help=(
            "Whether to filter to high coverage intervals for the sex ploidy and karyotype inference. "
            "Use only intervals that are considered high coverage across all platforms."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-platform-size",
        help=(
            "Required size of a platform to be considered in '--high-cov-by-platform-all'. Only platforms that "
            "have # of samples > 'min_platform_size' are used to determine intervals that have a high coverage across "
            "all platforms."
        ),
        type=int,
        default=100,
    )
    parser.add_argument(
        "--variant-depth-only-x-ploidy",
        help=(
            "Whether to use depth of variant data for the x ploidy estimation instead of the default behavior that "
            "will use reference blocks."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--variant-depth-only-y-ploidy",
        help=(
            "Whether to use depth of variant data for the y ploidy estimation instead of the default behavior that "
            "will use reference blocks."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--x-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the "
            "interval_coverage MT!"
        ),
        type=int,
        default=10,
    )
    parser.add_argument(
        "--y-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the "
            "interval_coverage MT!"
        ),
        type=int,
        default=5,
    )
    parser.add_argument(
        "--norm-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on the normalization autosome. This field must "
            "be in the interval_coverage MT!"
        ),
        type=int,
        default=20,
    )
    parser.add_argument(
        "--prop-samples-x",
        help="Proportion samples at specified coverage '--x-cov' to determine high coverage intervals on chromosome X.",
        type=float,
        default=0.80,
    )
    parser.add_argument(
        "--prop-samples-y",
        help="Proportion samples at specified coverage '--y-cov' to determine high coverage intervals on chromosome Y.",
        type=float,
        default=0.35,
    )
    parser.add_argument(
        "--prop-samples-norm",
        help=(
            "Proportion samples at specified coverage '--norm-cov' to determine high coverage intervals on the "
            "normalization chromosome specified by '--normalization-contig'."
        ),
        type=float,
        default=0.85,
    )
    parser.add_argument(
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'. Only "
            "used if '--test' is set and final coverage MT and/or platform assignment HT are not already written."
        ),
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    parser.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or 50 bp. Only used if '--test' is "
            "set and final coverage MT and/or platform assignment HT are not already written."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
