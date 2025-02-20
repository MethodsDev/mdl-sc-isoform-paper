# file of constants for consistent naming across the notebooks

# shortread PBMC samples
SHORTREAD_KEYS = ["pipseq_pbmc", "10x_3p_pbmc", "10x_5p_pbmc"]

# nice names for the shortread samples
SHORTREAD_NAMES = dict(zip(SHORTREAD_KEYS, ["PIPseq", "10x 3'", "10x 5'"]))


# consistent colors, used for both short and long read data
SAMPLE_COLORS = {"PIPseq": "#80BC42", "10x 3'": "#006DB6", "10x 5'": "#F36C3E"}

# monomer (long-read but not arrayed) samples
MONOMER_KEYS = {
    1: ("PIPseq", "0.8x"),
    2: ("PIPseq", "0.6x"),
    3: ("10x 3'", "purified"),
    4: ("10x 5'", "purified"),
    5: ("10x 3'", "unpurified"),
    6: ("10x 5'", "unpurified"),
}

# mas-seq (long-read, arrayed) samples
MASSEQ_KEYS = {
    1: ("PIPseq", "0.8x"),
    2: ("PIPseq", "0.6x"),
    3: ("10x 3'",),
    4: ("10x 5'",),
}


# terminal-friendly versions of the mas-seq sample names
MASSEQ_FILENAMES = {1: "pipseq_8x", 2: "pipseq_6x", 3: "10x_3p", 4: "10x_5p"}
