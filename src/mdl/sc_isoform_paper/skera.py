import csv
import itertools

import numpy as np


def read_ligation_csv(ligation_file, arraysize):
    counts = np.zeros((arraysize + 1, arraysize + 1), dtype=int)

    with open(ligation_file) as fh:
        for row in csv.DictReader(fh):
            counts[int(row["adapter_2"]), int(row["adapter_1"])] = int(row["ligations"])

    return counts


def read_length_csv(csv):
    used_zmws = set()
    len_to_concat = []
    mlen_to_concat = []

    with open(csv) as fh:
        for line in itertools.islice(fh, 1, None):
            zmw, hifi_rl, d_rl, concat = [int(x) for x in line.rstrip().split(",")]

            mlen_to_concat.append((d_rl, concat))
            if zmw in used_zmws:
                continue
            len_to_concat.append((hifi_rl, concat))
            used_zmws.add(zmw)

    return len_to_concat, mlen_to_concat
