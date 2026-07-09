import json
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# need to download human genome, see details in ./data/human_genome/
PATH_TO_GENOME = "./data/human_genome/GCF_000001405.25_GRCh37.p13_genomic.fna"
PATH_TO_JSON_OUT = "./data/triplet_counts_GRCh37.json"

NUCL_SET = set("ACGTacgt")


def is_appropriate_triplet(triplet: str) -> bool:
    return len(set(triplet).difference(NUCL_SET)) == 0


def main():
    triplet_counts = defaultdict(int)
    fasta = SeqIO.parse(PATH_TO_GENOME, "fasta")
    rec: SeqRecord = None
    for rec in fasta:
        print("Processing...", rec.description, file=sys.stderr)
        seq = str(rec.seq)
        # iterate over triplets with window=1
        for i in range(len(seq) - 2):
            triplet = seq[i: i + 3]
            if is_appropriate_triplet(triplet):
                triplet_counts[triplet] += 1

    print(triplet_counts, file=sys.stderr)
    with open(PATH_TO_JSON_OUT, "w") as fout:
        json.dump(triplet_counts, fout)


if __name__ == "__main__":
    main()
