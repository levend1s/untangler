
import argparse
import sys

from untangler_modules import __version__


def version():
    """Print the installed untangler version."""
    print(__version__)


def parse_bool(value):
    if isinstance(value, bool):
        return value
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes", "y", "on"}:
        return True
    if normalized in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(
        f"expected a boolean value, got {value!r}"
    )

from untangler_modules import cluster_transcripts


def build_parser():
    parser = argparse.ArgumentParser(
        prog="untangler",
        description="Untangler – cluster modBAM reads by transcript level features."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
        help="show the installed untangler version and exit",
    )

    subparsers = parser.add_subparsers(
        title="commands",
        dest="command"
    )

    # ---- cluster_transcripts command ----
    cluster_transcripts_parser = subparsers.add_parser("cluster_transcripts", help="cluster transcripts based on modification status")
    cluster_transcripts_parser.add_argument("-i", "--input",              required=False,  help="input file listing coverage data.")
    cluster_transcripts_parser.add_argument("-t", "--mod_prob_threshold", required=False, type=float, default=0.95, help="input file listing coverage data")
    cluster_transcripts_parser.add_argument("-s", "--strand",             required=False,  help="output filtered bam file.")
    cluster_transcripts_parser.add_argument("-m", "--mod",                required=False,  default="m6A", help="output filtered bam file.")
    cluster_transcripts_parser.add_argument("--exclude_tol",  required=False,  type=int, default=0, help="filter for only reads without m6A modification at these positions")
    cluster_transcripts_parser.add_argument("-a", "--annotation",         required=False,  help="annotation file (GFF)")
    cluster_transcripts_parser.add_argument("-g", "--genome",         required=False, help="genome for igv reports")
    cluster_transcripts_parser.add_argument("--ids", required=False, nargs="*", help="input file listing coverage data.")
    cluster_transcripts_parser.add_argument("--type", required=False, help="input file listing coverage data")

    cluster_transcripts_parser.add_argument("-p", "--padding", required=False, type=int, default=0, help="coverage padding")
    cluster_transcripts_parser.add_argument("-o", "--outfile", required=True, help="output directory")
    cluster_transcripts_parser.add_argument("-f", "--force", action="store_true", help="delete and recreate the output directory if it already exists")
    cluster_transcripts_parser.add_argument("-b", "--bamfile", required=False, help="input bam file with reads to cluster")
    cluster_transcripts_parser.add_argument("-d", "--min_deletion_length", required=False, type=int, default=30, help="min deletion length to consider as a variant indicator")
    cluster_transcripts_parser.add_argument("--cluster_cols", required=False, type=str, default=None, help="comma-separated list of columns to cluster on")
    cluster_transcripts_parser.add_argument("--min_cluster_percent", required=False, type=float, default=None, help="minimum percentage of reads to form a cluster (default: 0.01)")
    cluster_transcripts_parser.add_argument("--min_cluster_size_bulk", required=False, type=int, default=40, help="minimum number of reads to form a cluster (default: 100)")
    cluster_transcripts_parser.add_argument("--min_reads_to_process", required=False, type=int, default=40, help="minimum number of reads to process (default: 40)")
    cluster_transcripts_parser.add_argument("--distance_threshold", required=False, type=float, default=0.5, help="distance threshold for clustering (default: 0.1)")
    cluster_transcripts_parser.add_argument("--show_dendrogram", action="store_true", help="if provided, show dendrogram of clustering results")
    cluster_transcripts_parser.add_argument("--min_feature_freq", required=False, type=float, default=0.1, help="minimum frequency of a feature to be considered for clustering (default: 0.01)")
    cluster_transcripts_parser.add_argument("--min_gene_overlap", required=False, type=float, default=0, help="minimum fraction of read overlapping gene to consider for clustering (default: 0.5)")
    cluster_transcripts_parser.add_argument("--lift_threshold", required=False, type=float, default=1.1, help="minimum lift score to consider a feature as a predictor of intron presence (default: 1.5)")
    cluster_transcripts_parser.add_argument("--feature_distance_threshold", required=False, type=int, default=200, help="maximum distance between a feature and an intron to consider for clustering (default: 100)")
    cluster_transcripts_parser.add_argument("--hide_dendrogram_labels",  action="store_true", help="if provided, hide the labels on the dendrogram (default: False)")
    cluster_transcripts_parser.add_argument("--call_mods",  action="store_true", help="if provided, call modifications on the reads (default: False)")
    cluster_transcripts_parser.add_argument("--generate_igv_reports",  action="store_true", help="if provided, generate IGV reports for each cluster (default: False)")
    cluster_transcripts_parser.add_argument("--distance_metric", required=False, type=str, default="jaccard", help="distance metric to use for clustering (default: jaccard)")
    cluster_transcripts_parser.add_argument("--linkage_method", required=False, type=str, default="average", help="linkage method to use for clustering (default: average)")
    cluster_transcripts_parser.add_argument("--predict", required=False, type=str, default=None, help="column to use for prediction (default: None)")
    cluster_transcripts_parser.add_argument("--show_alignments", required=False, type=parse_bool, default=True, help="show alignment of reads in each cluster (True/False; default: True)")

    cluster_transcripts_parser.set_defaults(func=cluster_transcripts.cluster_transcripts)


    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()

    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)

    args.func(args)

if __name__ == "__main__":
    main()


