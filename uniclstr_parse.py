#!/usr/bin/env python3
import argparse
from collections import defaultdict

def parse_faa(faa_file):
    seqs = {}
    lengths = {}
    current = None
    with open(faa_file) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                seqs[current] = ""
            else:
                seqs[current] += line.strip()
    for sid, seq in seqs.items():
        lengths[sid] = len(seq)
    return seqs, lengths


def parse_clstr(clstr_file):
    reps = {}
    members = defaultdict(list)
    current_cluster = None

    with open(clstr_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[1])
            else:
                start = line.find(">") + 1
                end = line.find("...")
                pid = line[start:end]
                members[current_cluster].append(pid)
                if line.endswith("*"):
                    reps[current_cluster] = pid
    return reps, members


def extract_species(pid):
    return pid.split("|")[0]


def rename_headers_rep_only(seqs, reps):
    renamed = {}
    for cl, rep in reps.items():
        if rep in seqs:
            new_header = f"Uniclstr_{cl} {rep}"
            renamed[new_header] = seqs[rep]
    return renamed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
                        help="Input prefix (prefix.faa + prefix.clstr)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix")
    args = parser.parse_args()

    faa_file = args.input + ".faa"
    clstr_file = args.input + ".clstr"

    seqs, lengths = parse_faa(faa_file)
    reps, members = parse_clstr(clstr_file)

    tsv_out = args.output + ".tsv"
    with open(tsv_out, "w") as out:
        out.write("Uniclstr\tRepresentative\tRep_length\tSpecies\n")
        for cl in sorted(reps.keys()):
            rep = reps[cl]
            rep_len = lengths.get(rep, 0)
            species_set = { extract_species(pid) for pid in members[cl] }
            species_all = ";".join(sorted(species_set))
            out.write(f"Uniclstr_{cl}\t{rep}\t{rep_len}\t{species_all}\n")

    fasta_out = args.output + ".renamed.faa"
    renamed = rename_headers_rep_only(seqs, reps)
    with open(fasta_out, "w") as out:
        for header, seq in renamed.items():
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")


if __name__ == "__main__":
    main()
