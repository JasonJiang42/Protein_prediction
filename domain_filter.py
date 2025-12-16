#!/usr/bin/env python3
import argparse
import subprocess
from collections import defaultdict

def run_hmmscan(fasta, out, pfam):
    subprocess.run([
        "hmmscan",
        "--domtblout", out,
        pfam,
        fasta
    ], check=True)

def parse_domtbl(domtbl):
    doms = defaultdict(list)
    with open(domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            doms[fields[0]].append(fields[3])
    return doms

def arch(domlist):
    return ":".join(domlist)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", required=True, help="query FASTA")
    parser.add_argument("-i", "--input", required=True, help="candidate proteins FASTA")
    parser.add_argument("-p", "--pfam", required=True, help="path to Pfam-A.hmm")
    parser.add_argument("-o", "--output", required=True, help="output matched list")
    args = parser.parse_args()

    q_domtbl = "query.domtbl"
    c_domtbl = "candidates.domtbl"

    run_hmmscan(args.query, q_domtbl, args.pfam)
    run_hmmscan(args.input, c_domtbl, args.pfam)

    qdom = parse_domtbl(q_domtbl)
    cdom = parse_domtbl(c_domtbl)

    qid = list(qdom.keys())[0]
    qarch = arch(qdom[qid])

    with open(args.output, "w") as out:
        for pid, domlist in cdom.items():
            if arch(domlist) == qarch:
                out.write(f"{pid}\t{qarch}\n")

if __name__ == "__main__":
    main()