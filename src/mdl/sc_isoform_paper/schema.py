import os
import fastjsonschema


_path_validator = {"an existing path": os.path.exists}


def _path_exists(description: str = None):
    d = dict(type="string", format="an existing path")
    if description is not None:
        d["description"] = description
    return d


def _nonneg_int(description: str, default: int = None):
    d = dict(description=description, type="integer", minimum=0)
    if default is not None:
        d["default"] = default
    return d


def _pos_int(description: str, default: int = None):
    d = dict(description=description, type="integer", exclusiveMinimum=0)
    if default is not None:
        d["default"] = default
    return d


BARCODE_SCHEMA = {
    "type": "object",
    "properties": {
        "sample_type": {
            "description": "the scRNAseq library format",
            "enum": ["PIPseq", "10x 3'", "10x 5'"],
        },
        "barcode_file": _path_exists("gzipped text file containing barcode whitelist"),
        "umi_size": _pos_int("size of the UMI sequence"),
        "buffer_size": _pos_int("# bases to search for barcode and UMI"),
        "bam_paths": {
            "description": "list of marti-annotated BAM files to process",
            "type": "array",
            "items": _path_exists(),
        },
    },
    "required": ["sample_type", "barcode_file", "bam_paths"],
    "additionalProperties": False,
}


_PRIMING_PARAM_SCHEMA = {
    "type": "object",
    "properties": {
        "feature_pre": _nonneg_int("upstream bases for feature overlap", 5),
        "feature_post": _nonneg_int("downstream bases for feature overlap", 5),
        "motif_pre": _nonneg_int("upstream bases for polyA motifs", 20),
        "motif_post": _nonneg_int("downstream bases for polyA motifs", 30),
        "pas_pre": _nonneg_int("upstream bases for annotated polyA sites", 5),
        "pas_post": _nonneg_int("downstream bases for annotated polyA sites", 20),
        "polya_window": _nonneg_int("downstream bases for genomic polyA calls", 20),
        "polya_max_n": _nonneg_int("max As allowed in polyA window", 12),
        "polya_max_len": _nonneg_int("max consecutive As allowed in polyA window", 6),
    },
    "additionalProperties": False,
}

PRIMING_SCHEMA = {
    "type": "object",
    "properties": {
        "reference_fasta": _path_exists("fasta file for reference genome"),
        "reference_gtf": _path_exists("reference GTF, optionally converted to parquet"),
        "polya_motif_file": _path_exists("text file containing polyA motifs"),
        "polya_annotations": _path_exists("reference GTF of annotated polyA sites"),
        "priming_parameters": _PRIMING_PARAM_SCHEMA,
        "bam_paths": {
            "description": "list of BAM files aligned with minimap2",
            "type": "array",
            "items": _path_exists(),
        },
        "isoquant_paths": {
            "description": "IsoQuant output folder for each BAM",
            "type": "array",
            "items": _path_exists(),
        },
        "output_tag": {
            "description": "priming tag(s) to add to output",
            "type": "array",
            "items": {
                "description": "full classification and/or boolean tag",
                "enum": ["full", "simple"],
            },
            "default": ["full"],
            "minItems": 1,
            "uniqueItems": True,
        },
    },
    "required": [
        "reference_fasta",
        "reference_gtf",
        "polya_motif_file",
        "polya_annotations",
        "priming_parameters",
        "bam_paths",
    ],
    "additionalProperties": False,
}

barcode_validator = fastjsonschema.compile(BARCODE_SCHEMA, formats=_path_validator)
priming_validator = fastjsonschema.compile(PRIMING_SCHEMA, formats=_path_validator)
