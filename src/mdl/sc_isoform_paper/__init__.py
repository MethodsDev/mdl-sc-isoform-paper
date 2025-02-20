from datetime import date

# IUPAC basepair codes
_BASEPAIRS = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
        ".": ".",
        "-": "-",
    }
)


def rev_comp(seq: str):
    """reverse complement a sequence"""
    return seq.translate(_BASEPAIRS)[::-1]


class _Today:
    """object that returns today's date when formatted in a string"""

    def __str__(self):
        return date.today().strftime("%Y%m%d")


today = _Today()
