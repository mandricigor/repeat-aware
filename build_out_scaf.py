
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import shutil
from fastaq import sequences, utils, intervals, tasks



prefix = sys.argv[3]
similarity_level = sys.argv[4]
smallest_contig_length = sys.argv[5]

utils.syscall(' '.join(['nucmer', '-p', prefix, sys.argv[1], sys.argv[2], '--maxmatch']))
utils.syscall(' '.join(['delta-filter', '-i %s -l %s -r' % (similarity_level, smallest_contig_length), prefix + ".delta", '>', prefix + ".filter"]))
utils.syscall(' '.join(['show-coords', '-dTlro', prefix + ".filter", '>', prefix + ".coords"]))


output_scaffolds = {}
for record in SeqIO.parse(sys.argv[1], "fasta"):
    output_scaffolds[record.id] = str(record.seq)

ref_contigs = {}
for record in SeqIO.parse(sys.argv[2], "fasta"):
    ref_contigs[record.id] = str(record.seq)


#prefix = "prefix"



with open(prefix + ".coords") as f:
    a = f.readlines()

a = a[4:] # drop the header lines from nucmer hit file
a = list(map(lambda x: x.strip().split(), a))


contig_hits = [x[12] for x in a]
contig_hits_counts = Counter(contig_hits)

# groups is the list of groups by reference
groups = [[a[0]]]
linegroups = [[0]]
for i in range(1, len(a)):
    if a[i][11] == groups[-1][-1][11]:
        groups[-1].append(a[i])
        linegroups[-1].append(i)
    else:
        groups.append([a[i]])
        linegroups.append([i])


# because we may have that some contigs 
allgroups = []
alllinegroups = []
for linegroup, group in zip(linegroups, groups):
    group2 = [[group[0]]]
    linegroup2 = [[linegroup[0]]]
    for j in range(1, len(group)):
        if group[j][12] == group2[-1][-1][12]:
            group2[-1].append(group[j])
            linegroup2[-1].append(linegroup[j])
        else:
            group2.append([group[j]])
            linegroup2.append([linegroup[j]])
    allgroups.append(group2)
    alllinegroups.append(linegroup2)


contig_map = {}
used_contigs = set()

count = 0
scaffolds = []
chosen_lines = []
for alllinegroup, group in zip(alllinegroups, allgroups):
    scaf = []
    chosen = []
    for linegroup2, group2 in zip(alllinegroup, group):
        if len(group2) == 1 and len(group2[0]) == 13:
            # if the contig does not match completely -> skip it
            # we do not need it anymore
            name = group2[0][12]
            if int(group2[0][1]) - int(group2[0][0]) > 0.8 * int(group2[0][8]):
                scaf.append(name + ":::" + group2[0][10] + ":::" + str(int(group2[0][1]) - int(group2[0][0])))
                chosen.append(linegroup2)
                used_contigs.add(group2[0][12])
            #print(aaa)
            #print ("CONTINUE")
            #continue
        elif len(group2) == 1 and len(group2[0]) == 14 and group2[0][-1] != "[CONTAINED]":
            # it has probably a full match
                name = group2[0][12]
                scaf.append(name + ":::" + group2[0][10] + ":::" + str(int(group2[0][1]) - int(group2[0][0])))
                chosen.append(linegroup2)
                used_contigs.add(group2[0][12])
        else:
            # here we may have contigs corresponding to circular genomes
            # if it has both [BEGINS] and [ENDS] -> we just split it
            # or we still have mis-assembly - just merge all hits
            if int(group2[-1][1]) - int(group2[0][0]) > 0.8 * int(group2[0][8]):
                name = group2[0][12]
                scaf.append(name + ":::" + group2[0][10] + ":::" + str(int(group2[-1][1]) - int(group2[0][0])))
                chosen.append(linegroup2)
                used_contigs.add(group2[0][12])
            else:
                # check if this is the only hit
                #print ("DECI SUKA")
                name = group2[0][12]
                if group2[0][12] in contig_hits_counts and contig_hits_counts[group2[0][12]] == 1:
                    scaf.append(name + ":::" + group2[0][10] + ":::" + str(int(group2[-1][1]) - int(group2[0][0])))
                    used_contigs.add(group2[0][12])
                #chosen.append(linegroup2)
        count += 1
    if scaf:
        chosen_lines.append(chosen)
        scaffolds.append(scaf)


#print("COUNT:" + str(count))

unused_contigs = set(ref_contigs.keys()) - used_contigs


for refcont in unused_contigs:
    for outscaf, outseq in output_scaffolds.items():
        if outseq in ref_contigs[refcont]:
            #print(refcont + " " +  outscaf)
            break
        elif outseq in str(Seq(ref_contigs[refcont]).reverse_complement()):
            #print (refcont + " " + outscaf)
            break


distances = []
for cluster in chosen_lines:
    dist = [0]
    for group_number in range(len(cluster) - 2, -1, -1):
        dist.insert(0, int(a[cluster[group_number + 1][0]][0]) - int(a[cluster[group_number][-1]][1]))
    distances.append(dist)


rev_contig_map = {}
for x, y in contig_map.items():
    rev_contig_map[y] = x


with open(prefix + ".scaf", "w") as f:
    counter = 1
    for dist, scaffold in zip(distances, scaffolds):
        f.write(">scaffold_%s\n" % counter)
        for d, scaf in zip(dist, scaffold):
            f.write("%s:::%s\n" % (scaf, d))
        counter += 1
    for contig in unused_contigs:
        f.write(">scaffold_%s\n" % counter)
        f.write("%s:::1:::%s:::0\n" % (contig, len(ref_contigs[contig])))
        counter += 1









