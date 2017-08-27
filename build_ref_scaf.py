#!/usr/bin/env python3

import argparse
import os
import copy
import shutil
import sys
from fastaq import sequences, utils, intervals, tasks


from Bio import SeqIO
import sys



# check required nucmer programs are in path
progs = ['nucmer', 'delta-filter', 'show-coords']
not_in_path = [p for p in progs if shutil.which(p) is None]
    
if len(not_in_path):
    print('Error! Need these programs to be in your path:', file=sys.stderr)
    print('\n'.join(not_in_path), file=sys.stderr)
    print('Cannot continue', file=sys.stderr)
    sys.exit(1)


def nucmer_file_reader(fname):
    f = utils.open_file_read(fname)
    in_header = True

    for line in f:
        if in_header:
            if line.startswith('['):
                in_header = False
            continue
        yield NucmerHit(line)

    utils.close(f)


class NucmerHit:
    def __init__(self, line):
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  MAL1    NODE_25757_length_24482_cov_18.920391   [CONTAINS]
        try:
            l = line.rstrip().split('\t')
            self.ref_start = int(l[0])
            self.ref_end = int(l[1])
            self.qry_start = int(l[2])
            self.qry_end = int(l[3])
            self.hit_length_ref = int(l[4])
            self.hit_length_qry = int(l[5])
            self.percent_identity = float(l[6])
            self.ref_length = int(l[7])
            self.qry_length = int(l[8])
            self.frame = int(l[9])
            self.strand = int(l[10])
            self.ref_name = l[11]
            self.qry_name = l[12]

            if len(l) == 14:
                self.tag = l[13][1:-1]
            else:
                self.tag = None
        except:
            print('Error reading this nucmer line:\n' + line, file=sys.stderr)

    def __str__(self):
         return str(self.qry_name) + " : " + str(self.ref_start) + " - " + str(self.ref_end)



parser = argparse.ArgumentParser(
    description = 'Takes contigs and a reference sequence. Makes a new fasta file of the contigs, but they are now perfect sequences by using the reference instead',
    usage = '%(prog)s [options] <contigs.fa> <reference.fa> <outprefix> <smallest.contig.length> <similarity.level>')
parser.add_argument('--nucmer_options', help='Options when running nucmer [%(default)s]', default='')
parser.add_argument('contigs_fa', help='Name of contigs fasta file', metavar='contigs.fa')
parser.add_argument('ref_fa', help='Name of reference fasta file', metavar='reference.fa')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('min_seq_length', type=int, help='Minimum length of contig to output', metavar='smallest.contig.length')
parser.add_argument('sim_level', type=int, help='Minimum similarity level of nucmer hits to consider', metavar='similarity.level')
options = parser.parse_args()

ref_seqs = {}
tasks.file_to_dict(options.ref_fa, ref_seqs)

nucmer_out_prefix = options.outprefix
nucmer_out_delta = nucmer_out_prefix + '.delta'
nucmer_out_filter = nucmer_out_prefix + '.filter'
nucmer_out_coords = nucmer_out_prefix + '.coords'
smallest_contig_length = options.min_seq_length
similarity_level = options.sim_level

# load contigs into memory
sw = {}
handle = open(options.contigs_fa, "rU")
for record in SeqIO.parse(handle, "fasta"):
    sw[record.id] = record.seq.tostring()
handle.close()





# run nucmer of contigs vs ref
utils.syscall(' '.join(['nucmer', options.nucmer_options, '-p', nucmer_out_prefix, options.ref_fa, options.contigs_fa, '--maxmatch']))
utils.syscall(' '.join(['delta-filter', '-i %s -l %s -r' % (similarity_level, smallest_contig_length), nucmer_out_delta, '>', nucmer_out_filter]))
utils.syscall(' '.join(['show-coords', '-dTlro', nucmer_out_filter, '>', nucmer_out_coords]))

# load hits into hash. key=ref_name, value=another hash with key=qry_name, value=list of hit positions in that ref seq
nucmer_hits = {}
contigs_to_print = {}




nucmer_reader = nucmer_file_reader(nucmer_out_coords)

# this dictionary contains hits for each contig per reference
hit_ref_dict = {}

for hit in nucmer_reader:
    hit_dict = hit_ref_dict.setdefault(hit.qry_name, {})
    ref_list = hit_dict.setdefault(hit.ref_name, [])
    ref_list.append(hit)



qry_segments = {}
for qry_name, qry_dict in hit_ref_dict.items():
    # determine if we have a full match, if so -> just add it to the list of contigs
    contains = False
    good_match = None
    for ref_name, ref_list in qry_dict.items():
        contains_contigs = [hit for hit in ref_list if hit.tag == "CONTAINS"]
        if contains_contigs:
            contains = True
            good_match = contains_contigs[0]
            break
    if contains:
        # if we have a hit for a contig which is completely contained 
        # then we just consider it -> we want to take biggest hit
        # and ignore the rest of them
        a, b = (good_match.qry_start, good_match.qry_end)
        if a > b:
            qry_segments[qry_name] = [(b, a)]
        else:
            qry_segments[qry_name] = [(a, b)]
    else:
        segments = []
        for ref_name, ref_list in qry_dict.items():
            hits = sorted(ref_list, key=lambda z: z.ref_start)
            # merge all overlapping hits for each reference
            merging = [[0]]
            for i in range(1, len(hits)):
                if hits[i].ref_start <= hits[merging[-1][-1]].ref_end + 1:
                    if hits[i].ref_end > hits[merging[-1][-1]].ref_end:
                        merging[-1].append(i)
                else:
                    merging.append([i])
            for merge_group in merging:
                a, b = (hits[merge_group[0]].qry_start, hits[merge_group[-1]].qry_end)
                if a > b:
                    segments.append((b, a))
                else:
                    segments.append((a, b))
        qry_segments[qry_name] = segments



###################################################

preliminary = options.outprefix + '.fa'
# print the preliminary perfect contigs
with open(preliminary, "w") as f:
    for x, y in qry_segments.items():
        if len(y) > 1:
            counter = 1
            for u, v in y:
                f.write(">%s\n" % (x + "_" + str(counter)))
                f.write("%s\n" % sw[x][(u - 1):(v - 1)])
                counter += 1
        else:
            f.write(">%s\n" % x)
            f.write("%s\n" % sw[x][(y[0][0] - 1):(y[0][1] - 1)])





utils.syscall(' '.join(['nucmer', options.nucmer_options, '-p', nucmer_out_prefix, options.ref_fa, preliminary, '--maxmatch']))
utils.syscall(' '.join(['delta-filter', '-i %s -l %s -r' % (similarity_level, smallest_contig_length), nucmer_out_delta, '>', nucmer_out_filter]))
utils.syscall(' '.join(['show-coords', '-dTlro', nucmer_out_filter, '>', nucmer_out_coords]))











# load into memory the contigs
sw = {}
handle = open(preliminary, "rU")
for record in SeqIO.parse(handle, "fasta"):
    sw[record.id] = record.seq.tostring()
handle.close()


prefix = options.outprefix

# load coordinates
with open(nucmer_out_coords) as f:
    a = f.readlines()

a = a[4:] # drop the header lines from nucmer hit file
a = list(map(lambda x: x.strip().split(), a))

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
global_counter = 0

scaffolds = []
chosen_lines = []
for alllinegroup, group in zip(alllinegroups, allgroups):
    scaf = []
    chosen = []
    for linegroup2, group2 in zip(alllinegroup, group):
        if len(group2) == 1 and len(group2[0]) == 12:
            # if the contig does not match completely -> skip it
            # we do not need it anymore
            continue
        elif len(group2) == 1 and len(group2[0]) == 13:
            # it has probably a full match
            if group2[0][12] == "[CONTAINS]":
                if group2[0][12] in contig_map:
                    name = contig_map[group2[0][12]]
                else:
                    global_counter += 1
                    name = "contig_%s" % global_counter
                    contig_map[group2[0][12]] = name
                scaf.append(name + ":::" + group2[0][10] + ":::" + str(int(group2[0][1]) - int(group2[0][0])))
                chosen.append(linegroup2)
        else:
            # here we may have contigs corresponding to circular genomes
            # if it has both [BEGINS] and [ENDS] -> we just split it
            # or we still have mis-assembly - just merge all hits
            if int(group2[-1][1]) - int(group2[0][0]) > 0.97 * int(group2[0][8]):
                #print group2, "MMMMMMMMMMMM"
                if group2[0][12] in contig_map:
                    name = contig_map[group2[0][12]]
                else:
                    global_counter += 1
                    name = "contig_%s" % global_counter
                    contig_map[group2[0][12]] = name
                scaf.append(name + ":::" + group2[0][10] + ":::" + str(int(group2[-1][1]) - int(group2[0][0])))
                chosen.append(linegroup2)
    chosen_lines.append(chosen)
    scaffolds.append(scaf)


distances = []
for cluster in chosen_lines:
    dist = [0]
    for group_number in range(len(cluster) - 2, -1, -1):
        dist.insert(0, int(a[cluster[group_number + 1][0]][0]) - int(a[cluster[group_number][-1]][1]))
    distances.append(dist)


rev_contig_map = {}
for x, y in contig_map.items():
    rev_contig_map[y] = x



with open("%s.scaf" % prefix, "w") as f:
    counter = 1
    for dist, scaffold in zip(distances, scaffolds):
        f.write(">scaffold_%s\n" % counter)
        for d, scaf in zip(dist, scaffold):
            f.write("%s:::%s\n" % (scaf, d))
        counter += 1


with open("%s.fa" % prefix, "w") as f:
    all_contigs = set()
    for scaffold in scaffolds:
        for contig in scaffold:
            all_contigs.add(contig.split(":::")[0])
    scontigs = []
    for sc in all_contigs:
        scontigs.append((sc, int(sc.split("_")[1])))
    for contig, order in sorted(scontigs, key=lambda z: z[1]):
        f.write(">%s\n" % contig)
        f.write("%s\n" % sw[rev_contig_map[contig]])


