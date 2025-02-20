import itertools
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import click
import yaml

from mdl.log import verbosity_config_option
from bouncer import BarcodeSet

from . import barcode_extraction as bc_ext
from . import priming
from . import schema

log = logging.getLogger(__package__)


def _share_barcode_set(b_set):
    bc_ext.barcode_set = b_set


def schema_option(schema_dict):
    def print_schema(ctx: click.Context, param: str, value: str):
        if not value or ctx.resilient_parsing:
            return

        click.echo(yaml.dump(schema_dict, sort_keys=False))
        ctx.exit()

    return click.option(
        "--schema",
        help="Print the schema for the config file and exit.",
        is_flag=True,
        callback=print_schema,
        expose_value=False,
        is_eager=True,
    )


@click.command
@click.option(
    "--config-file", required=True, type=click.Path(exists=True, path_type=Path)
)
@click.option(
    "--processes",
    "-p",
    type=int,
    default=1,
    help="Number of processes to use.",
)
@click.option(
    "--output-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Where to write output BAMs. Default: next to input",
)
@verbosity_config_option(log, pkg=__package__)
@click.version_option()
@schema_option(schema.BARCODE_SCHEMA)
def extract_barcodes(config_file: Path, processes: int, output_dir: Path = None):
    """Extracts barcodes and UMIs from long reads of single-cell RNAseq data and adds
    them to the BAM tags. Requires BAMs tagged with marti annoations"""

    log.info("Loading config from %s", config_file)
    with open(config_file, "rt") as fh:
        config = schema.barcode_validator(yaml.safe_load(fh))

    barcode_set = BarcodeSet.load_from(Path(config["barcode_file"]))

    if config["sample_type"] == "PIPseq":
        tag_func = bc_ext.extract_3p_bc_umi
        buffer_size = 56
        umi_size = 12
    elif config["sample_type"] == "10x 3'":
        tag_func = bc_ext.extract_3p_bc_umi
        buffer_size = 29
        umi_size = 12
    elif config["sample_type"] == "10x 5'":
        tag_func = bc_ext.extract_5p_bc_umi
        buffer_size = 27
        umi_size = 10
    else:
        raise ValueError(f"Invalid sample type {config['sample_type']}")

    bam_paths = list(map(Path, config["bam_paths"]))
    tagged_bams = [bp.with_suffix(".tagged.bam") for bp in bam_paths]
    if output_dir is not None:
        tagged_bams = [output_dir / tb.name for tb in tagged_bams]

    log.debug(f"Creating process pool of {processes} processes")
    with ProcessPoolExecutor(
        processes, initializer=_share_barcode_set, initargs=(barcode_set,)
    ) as exc:
        tagged_reads = dict(
            exc.map(
                bc_ext.write_barcode_bam,
                bam_paths,
                tagged_bams,
                itertools.repeat(tag_func),
                itertools.repeat(config.get("umi_size", umi_size)),
                itertools.repeat(config.get("buffer_size", buffer_size)),
            )
        )
    log.info(f"Finishing annotating {len(config['bam_paths'])} BAM files")
    log.info("Wrote %s reads", sum(v + 1 for v in tagged_reads.values()))


@click.command
@click.option(
    "--config-file", required=True, type=click.Path(exists=True, path_type=Path)
)
@click.option(
    "--processes",
    "-p",
    type=int,
    default=1,
    help="Number of processes to use.",
)
@click.option(
    "--output-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Where to write output BAMs. Default: next to input",
)
@click.option(
    "--filter-isoquant-by-bam",
    is_flag=True,
    help="Use the name of each BAM file to filter the IsoQuant files.",
)
@verbosity_config_option(log, pkg=__package__)
@click.version_option()
@schema_option(schema.PRIMING_SCHEMA)
@click.pass_context
def annotate_priming(
    ctx: click.Context,
    config_file: Path,
    processes: int,
    output_dir: Path = None,
    filter_isoquant_by_bam: bool = False,
):
    """Annotates mapped BAM files with priming classifications based on a set
    of configurable heuristics"""

    log.info("Loading config from %s", config_file)
    with open(config_file, "rt") as fh:
        config = schema.priming_validator(yaml.safe_load(fh))

    log.info("Loading reference files")
    reference = priming.ReferenceRanges(
        reference_fasta=Path(config["reference_fasta"]),
        reference_gtf=Path(config["reference_gtf"]),
        polya_motif_file=Path(config["polya_motif_file"]),
        polya_annotations=Path(config["polya_annotations"]),
    )
    priming_classifier = priming.PrimingClassifier(**config["priming_parameters"])

    bam_paths = list(map(Path, config["bam_paths"]))
    anno_bams = [bp.with_suffix(".annotated.bam") for bp in bam_paths]
    if output_dir is not None:
        anno_bams = [output_dir / ab.name for ab in anno_bams]

    if "isoquant_paths" in config:
        if (n_b := len(bam_paths)) != (n_iq := len(config["isoquant_paths"])):
            log.error("IsoQuant paths are not one-to-one with BAM files!")
            log.error(f"Found {n_b} BAMs and {n_iq} IsoQuant paths")
            ctx.exit(1)
        isoquant_paths = map(Path, config["isoquant_paths"])
    else:
        isoquant_paths = itertools.repeat(None)

    log.debug(f"Creating process pool of {processes} processes")
    with ProcessPoolExecutor(processes) as exc:
        list(
            exc.map(
                priming_classifier.tag_bam_with_read_stats,
                bam_paths,
                anno_bams,
                itertools.repeat(reference),
                isoquant_paths,
                itertools.repeat(filter_isoquant_by_bam),
                itertools.repeat("full" in config["output_tag"]),
                itertools.repeat("simple" in config["output_tag"]),
            )
        )
    log.info(f"Finishing annotating {len(config['bam_paths'])} BAM files")
