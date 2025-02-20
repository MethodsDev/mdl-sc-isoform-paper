import itertools


# four tiers of 96 sequences used by Fluent PIPseq to create barcode whitelist

TIER_1 = """
AGAAACCA CCTTTACA TACCTCCC ACCTTCCC AAACTACA ATTACCTT GATTTCCC CTCCTCCA
TCTAAACT AAGTCCAA TCCGACAC GAGAAACC ACCCTCAA AGACCTCA GATTACTT CCACCTCT
TAACTTCT TTCCCTAT CTGTTTCC TCCTATAT CACTAACC CCCTGTTT TTGACCCA TCTATTCC
GTGTCACC AGGACACA CTTTGGAC CCTATTTA TAGTCTCT CTTTCACT ATCCCACC ATACTCTC
CAATTCTC ATTTCCAT CAAGGGTT GTCTTCCT CTGGGTAT CAAACATT CAGGTTGC GTCCTTGC
GATTGGGA TTGGGTCC TACCCTGC GAGGGTCA AGAGGTGC CTGTGACC GTCCACTA CTTAGTGT
GAGTGTAC GTTGTCCG TCTTTGAC CCTTTGTC GTGAACTC AAGGGACC AATACATC AAACAAAC
ACTACCCG GATGTGGC ACCCATGC ACCAACCC AAAGAGGC AAGTTGTC GAATCCCA CTTTATCC
GTAAACAA CTTCTACG CCATCCAC CCTCATGA CTAGACTA ACCAGTTT AGTTGAAC AGTTTGTA
CCCTCTTG GAGGAGTG TCCCTGGA AAATTCCG GAAATACG AGTCACAA TACTGAAT AGACGAGG
CCTACGCT AAACCGCC AATATGAC GACACCTG CTGTTGTG CTAACGCC TTCACTGG GTCTAATC
TATGTGAA GTGAGGCA TATCTGTC GTGGTGCT ATACACCC CCCTTGCA TATCCACG GCCTGGTA""".split()

TIER_2 = """
AGGAAA AGGTAA AGTGGA ATGTTG GGTTTC GTAGAG GTTAGT GTTTGG
TAGCGA TATTGG TGGGTT TTGGTA AAGAGA AAAGTG TAAGGC AAAGGC
AATAGC TAAGCC TATGCC GTTGCT CAGTTG GAAAGG CACAAG TACAGA
GAAGAA ACAAAG AGAAGG GCGTTT TGAAAG TGAGAA GATGAA CACGAA
ACGGTT GATTTC TAGTCT TTTCTC CGCAAA CATCTA GGTCTA AAGGTG
AAAGAC AGAAAC TTAACG TGAACC AGTTAC AAACCG TAACCC GCACTA
AGACGT ACATGT ATCAAC TTCGAA GACAAT TGCTTT TGCTAG GTGATC
GATATG GTAATC GAAATC GCTGTA TGTATC CTGAAG ATGCAC GTACAA
AAACAC AAGCAC GAACAG GTTCAC ACCTTT AACTGA CCGTAT ATCTGA
TCAAAG TCGATT GCTAAG GAGATA CTGGTA CTTGTT CTTTAG CTTCGA
CTAAAG CTATGG TCTGTG CACATT TCAGTC CCAAAT AATACC ACTTCC
ACAACC CCTAAT ACAGCA CTTGAC CAATAC GCTCTT TTGGCA TCTACC""".split()

TIER_3 = """
AAGGTG AGGAAA AGGTAA AGTGGA ATGCAC ATGTTG GGTTTC GTAATC
GTACAA GTAGAG GTGATC GTTAGT GTTTGG TAGAAC TAGCGA TATTGG
TGGGTT TGTAAC TGTATC TTAACG TTGGCA TTGGTA AAGAGA AAAGTG
TAAGGC AAACCG AAACAC AAAGAC AAAGGC AATAGC AAGCAC AGTTAC
TAAGCC TATGCC CTGAAG CTGGTA GTTGCT GCTGTA CAGTTG CTTGAC
CTTGTT CTTTAG CTAAAG CTATGG GCTAAG TAACCC AATACC ACAGCA
ACATGT ACTTCC TCAGTC CAAACA CAATAC TGCTTT TGCTAG TCTGTG
CTTCCA CCTAAT AACTGA GAAAGG GAACAG CACAAG CACATT TACAGA
ACCATA TCCATA GTTCAC TCTACC TTCCAG GAAGAA ACAAAG ACAACC
CCAAAT GACAAT GCACTA ATCAAC GAAACC GAAATC GATATG AGAAAC
AGAAGG AGACGT GATACC TCAAAG ACCTTT AAGCTC GCGTTT CTTCGA
ATCTGA TGAAAG TGAACC TGAGAA ATCTTC AAACTC TACTCA CACGAA""".split()

TIER_4 = """
CCTATTTA CAGTTTAA ATACACCC ATCCCACC TCCTATAT TTAGCAAT ACCAACCC GCCAACAT
ACCAGTTT AGTAGTTA CCTTTACA ACAAGTAG ACACCAAG ACCCATGC GGCCCAAT GGCTATAA
AGCATGCC CACTAACC ATGCATAT TTGCATTC AATAAGGA TTCTAGGA AGTAATGG ACTAATTG
ATCAGGGA CTCCCAAA GACACAAA CCATATGA ATATGCAA GCTAAGTT TGCACCAG CAAACATT
AAACTACA ATCCTAGT TACCTAAG TGGCTAGT CACAACCT ATAACAGG GGCAAGGT GGTTAGGG
GTCAGGTT CTCAAACA AATATGAC AACAGAAC ACCACAGA ACTAGAGC GATGCAGA GTCAAGAG
TCCAGAAG GGTTACAC AAACTGTG TAAGGGCC TAGTAGCC ACAGGCCA CAAGGAAT CAAGGTAC
GCTATGGG AACAAATG CCTTTGTC CCCTGTTT TTTGCCAG GCCTGGTA GTCTGGAA TCTGATTT
ACAAAGAT AACTGCCT CTCACATC GTCTAATC AAATAGCA AATGTATG CGTGTACA GATGTGGC
TATGTGAA CGTGGGAT GATGGTTA AAAGAGGC ACAGATAA ATAGATGT CCAGACAG GAAGATAT
TGGGAATT AGTTTGTA GGAGGTTT GGAGTAAG ACCTGAAG TACTGAAT CAAGGGTT CTGTTAAA
CTGTTTCC GAGTGTAC TAATGTGG CTGGGTAT TTGGGTCC ACATGGAC ATATGGGT GAAAGACA""".split()


# map from barcode to index, used for PIPseeker conversion
TIER_1_DICT = {s: i for i, s in enumerate(TIER_1)}
TIER_2_DICT = {s: i for i, s in enumerate(TIER_2)}
TIER_3_DICT = {s: i for i, s in enumerate(TIER_3)}
TIER_4_DICT = {s: i for i, s in enumerate(TIER_4)}

# map from base to binary encoding and vice versa
BASE_CODE = {"A": "00", "C": "01", "G": "10", "T": "11"}
BIN_CODE = {("0", "0"): "A", ("0", "1"): "C", ("1", "0"): "G", ("1", "1"): "T"}


# NOTE: Over time, PIPseq changed how it encodes barcodes. In the first version, the raw
# sequence was encoded as four integers, each indexing into one of the whitelists, and
# these integers were individually encoded as groups of four bases. In the later version
# the whitelist is numbered from 0 to 96**4 - 1 and this number is encoded with 16 bases


def _bin_to_base(i: int) -> str:
    """Translate from a binary string to the bases used by pipseeker"""
    return "".join(
        BIN_CODE[b]
        for b in itertools.islice(
            itertools.pairwise(bin(i)[2:].zfill(8)), None, None, 2
        )
    )


def _base_to_bin(bc: str) -> tuple[int, int, int, int]:
    """Translate from the pipseeker barcode to the four indexes into the sequence lists"""
    return tuple(
        int("".join(BASE_CODE[c] for c in bc[i : i + 4]), base=2)
        for i in range(0, 16, 4)
    )


def _tuple_to_sequence(t1: int, t2: int, t3: int, t4: int) -> str:
    """Index into the four whitelists and combine with the linker sequences"""
    return f"{TIER_1[t1]}ATG{TIER_2[t2]}GAG{TIER_3[t3]}TCGAG{TIER_4[t4]}"


def barcode_to_int(bc: str) -> int:
    """Given an artificial barcode, converts it to an integer that combines the indexes
    into the four different whitelists. Used in PIPseeker v3.3 and useful for storage"""
    return int("".join(BASE_CODE[c] for c in bc), base=2)


def barcode_to_sequence(bc: str) -> str:
    """Given the artificial barcode from PIPseeker v3, return the real sequence"""
    return _tuple_to_sequence(*_base_to_bin(bc))


def barcode_to_sequence_v33(bc: str) -> str:
    """Given the artificial barcode from PIPseeker v3.3+, return the real sequence"""
    i = barcode_to_int(bc)
    t1 = i % 96
    i //= 96
    t2 = i % 96
    i //= 96
    t3 = i % 96
    i //= 96
    t4 = i % 96

    return _tuple_to_sequence(t1, t2, t3, t4)


def _sequence_to_tuple(seq: str) -> tuple[int, int, int, int]:
    """Given the real sequence, return the indices into the PIPseq whitelists"""
    t1 = TIER_1_DICT[seq[0:8]]
    t2 = TIER_2_DICT[seq[11:17]]
    t3 = TIER_3_DICT[seq[20:26]]
    t4 = TIER_4_DICT[seq[31:39]]

    return (t1, t2, t3, t4)


def sequence_to_barcode(seq: str) -> str:
    """Given the real sequence, return the PIPseeker v3 barcode"""
    return "".join(_bin_to_base(t) for t in _sequence_to_tuple(seq))


def sequence_to_barcode_v3(seq: str) -> str:
    """Given the real sequence, return the PIPseeker v3.3+ barcode"""
    i = sequence_to_int(seq)

    return "".join(
        BIN_CODE[b]
        for b in itertools.islice(
            itertools.pairwise(bin(i)[2:].zfill(32)), None, None, 2
        )
    )


def sequence_to_int(seq: str) -> int:
    """Given the real sequence, return an integer from 0 to 96**4 - 1"""
    t1, t2, t3, t4 = _sequence_to_tuple(seq)

    return 96 * (96 * (96 * t4 + t3) + t2) + t1


if __name__ == "__main__":
    # print out all the barcodes
    for t1, t2, t3, t4 in itertools.product(TIER_1, TIER_2, TIER_3, TIER_4):
        print(_tuple_to_sequence(t1, t2, t3, t4))
