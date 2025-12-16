#!/usr/bin/env python3
import argparse
import collections
import csv

# ------------------------------------------------------------
# Read GA thresholds from Pfam-A.hmm.dat
# ------------------------------------------------------------
def read_pfam_data(dat_file):
    Data = collections.namedtuple("Data", ["type","clan","ga_seq","ga_dom"])
    data = {}
    with open(dat_file) as fh:
        clan = None
        for line in fh:
            if line.startswith("#=GF ID"):
                hmm_name = line[10:].strip()
            elif line.startswith("#=GF TP"):
                typ = line[10:].strip()
            elif line.startswith("#=GF CL"):
                clan = line[10:].strip()
            elif line.startswith("#=GF GA"):
                vals = line[10:].strip().rstrip(";").split(";")
                ga_seq = float(vals[0])
                ga_dom = float(vals[1])
            elif line.startswith("//"):
                data[hmm_name] = Data(typ, clan, ga_seq, ga_dom)
                clan = None
    return data


# ------------------------------------------------------------
# Parse domtblout + GA filtering (fixed)
# ------------------------------------------------------------
def parse_domtbl_with_ga(domtbl, pfam_data):
    kept = []
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split()

            # Ensure domtblout has at least 23 fields
            if len(cols) < 23:
                continue

            # rebuild description field
            description = " ".join(cols[22:])
            cols = cols[:22] + [description]

            hmm = cols[0]
            if hmm not in pfam_data:
                continue

            try:
                score_seq = float(cols[7])
                score_dom = float(cols[13])
            except ValueError:
                continue

            ga = pfam_data[hmm]

            if score_dom >= ga.ga_dom and score_seq >= ga.ga_seq:
                kept.append(cols)

    return kept


# ------------------------------------------------------------
# Sort rows by ali_from (column 17) within each TARGET sequence (column 3)
# ------------------------------------------------------------
def sort_within_each_target(rows):
    groups = {}

    for r in rows:
        target = r[3]  # target (protein ID)
        groups.setdefault(target, []).append(r)

    # sort each target group
    for target in groups:
        groups[target].sort(key=lambda x: int(x[17]))  # ali_from ascending

    # rebuild final list
    out = []
    for target in groups:
        out.extend(groups[target])
    return out


# ------------------------------------------------------------
# Output CSV
# ------------------------------------------------------------
def write_csv(path, rows):
    header = [
        "target","target_acc","tlen","query","query_acc","qlen",
        "E_full","score_full","bias_full","domnum","domof",
        "c_E","i_E","score_dom","bias_dom",
        "hmm_from","hmm_to","ali_from","ali_to",
        "env_from","env_to","acc","description"
    ]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", required=True, help="input hmmscan results file (.domtbl)" )
    ap.add_argument("-d", required=True, help="Path to pfam.hmm.dat file for GA value")
    ap.add_argument("-o", required=True, help="output files")
    args = ap.parse_args()

    pfam = read_pfam_data(args.d)
    rows = parse_domtbl_with_ga(args.i, pfam)
    rows = sort_within_each_target(rows)
    write_csv(args.o, rows)


if __name__ == "__main__":
    main()