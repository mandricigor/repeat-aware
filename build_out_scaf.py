


from Bio import SeqIO
import re
import sys
from collections import Counter
import argparse



parser = argparse.ArgumentParser(description='Building output scaffolding')
parser.add_argument('--scafolds', dest="scaffolds", type=str, required=True, help="fasta file with output scaffolds")
parser.add_argument('--contigs', dest="contigs", type=str, required=True, help="fasta file with contigs")
parser.add_argument('--outscaf', dest="outscaf", type=str, required=True, help="output scaf file (output scaffolding)")

args = parser.parse_args()


contigs = {}
for record in SeqIO.parse(args.contigs, "fasta"):
    contigs[str(record.seq).upper()] = record.id + ":::1"
    contigs[str(record.seq.reverse_complement()).upper()] = record.id + ":::-1"

a = Counter(contigs.values())



scaffolds = []
for scaffold in SeqIO.parse(args.scaffolds, "fasta"):
    scaffolds.append(scaffold.seq)


scafs = []
for scaffold in scaffolds:
    myscaf = []
    metacontigs = re.split("[N|n]+", str(scaffold))
    distances = iter(map(len, re.findall("[N|n]+", str(scaffold))) + [0])
    for jj, meta in enumerate(metacontigs):
        #if jj % 1000 == 0:
        #    print jj
        #print jj, "MMM"
        contig = contigs.get(meta, "NA")
        if contig == "NA":
            contained = False
            good = None
            for key in contigs.keys():
                if meta.upper() in key:
                    contained = True
                    good = contigs[key]
                    break
            if contained:
                contig = good + ":::" + str(len(meta)) + ":::" + str(distances.next())
                myscaf.append(contig)
            else:
                # scaffold consists of several contigs sticked together
                intervals = {}
                for key in contigs.keys():
                    position = meta.upper().find(key)
                    if position != -1:
                        intervals[contigs[key] + ":::" + str(len(key))] = (position, position + len(key))
                intervals = sorted(intervals.items(), key=lambda x: x[1][0])
                remove = set()
                for i in range(len(intervals)):
                    for j in range(len(intervals)):
                        if i == j:
                            continue
                        int1 = intervals[i][1]
                        int2 = intervals[j][1]
                        if int1[0] >= int2[0] and int1[1] <= int2[1]:
                            # interval int1 is containing in another one - skip it
                            remove.add(i)
                intervals2 = []
                for i in range(len(intervals)):
                    if i not in remove:
                        intervals2.append(intervals[i])
                mymyscaf = []
                for i in range(len(intervals2) - 1):
                    mymyscaf.append(intervals2[i][0] + ":::" + str(intervals2[i + 1][1][0] - intervals2[i][1][1]))
                mymyscaf.append(intervals2[-1][0] + ":::" + str(distances.next()))
                myscaf.extend(mymyscaf)
        else:
            contig = contig + ":::" + str(len(meta)) + ":::" + str(distances.next())
            myscaf.append(contig)
    scafs.append(myscaf)


with open(args.outscaf, "w") as f:
    for i, scaf in enumerate(scafs):
        f.write(">scaffold_%s\n" % i)
        for c in scaf:
            f.write("%s\n" % c)



