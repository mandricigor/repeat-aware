

import networkx as nx
import sys
import subprocess
from Bio import SeqIO


def syscall(cmd):
    retcode = subprocess.call(cmd, shell=True)
    if retcode != 0:
        raise Exception("Error in system call. Command was:\n" + cmd)


ref_fa = sys.argv[1]
contigs_fa = sys.argv[2]
prefix = sys.argv[3]
similarity_level = int(sys.argv[4])
smallest_contig_length = int(sys.argv[5])


syscall(' '.join(['nucmer', '-p', prefix, ref_fa, contigs_fa, '--maxmatch', '--nosimplify']))
syscall(' '.join(['delta-filter', '-i %s -l %s ' % (similarity_level, smallest_contig_length), prefix + ".delta", '>', prefix + ".filter"]))
syscall(' '.join(['show-coords', '-dTlro', prefix + ".filter", '>', prefix + ".coords"]))

















sw = {}
for record in SeqIO.parse(contigs_fa, "fasta"):
    sw[record.id] = str(record.seq)


with open(prefix + ".coords") as f:
    a = f.readlines()


a = map(lambda x: x.strip().split(), a[4:])

# remove duplicate lines
aa = [a[0]]
for line in a[1:]:
    if line == aa[-1]:
        continue
    else:
        aa.append(line)

a = aa
# removed duplicate lines



# count the occurence of each contig
counts = {}
for line in a:
    if line[12] not in counts:
        counts[line[12]] = 0
    counts[line[12]] += 1
# found counts of each contig


refs = [[a[0]]]
for line in a[1:]:
    if line[11] == refs[-1][-1][11]:
        refs[-1].append(line)
    else:
        refs.append([line])




# find contigs contained inside others
contained = {}
junk_lines = []
for ref in refs:
    aa = [[ref[0]]]
    for line in ref[1:]:
        c1 = int(line[0])
        c2 = int(line[1])
        if c1 <= int(aa[-1][-1][1]) and c2 <= int(aa[-1][-1][1]):
            if line[12] not in contained:
                contained[line[12]] = 0
            contained[line[12]] += 1
            junk_lines.append(line)
        elif c1 <= int(aa[-1][-1][1]) and c2 > int(aa[-1][-1][1]):
            aa[-1].append(line)
        else:
            aa.append([line])



# try without filtering out
#for jl in junk_lines:
#    print jl
#
#a = [xxx for xxx in a if xxx not in junk_lines]












refs = [[a[0]]]
for line in a[1:]:
    if line[11] == refs[-1][-1][11]:
        refs[-1].append(line)
    else:
        refs.append([line])







refdict = {}
for ref in refs:
    refdict[ref[0][11]] = ref





contig_coords = {} # all coordinates for all contigs
coords = {}
for ref, lines in refdict.items():
    cdict = {}
    for i, line in enumerate(lines):
        c1 = int(line[0])
        c2 = int(line[1])
        if c1 not in cdict:
            cdict[c1] = []
        if c2 not in cdict:
            cdict[c2] = []
        cdict[c1].append((line[12] + ":::start:::" + str(i + 1), int(line[2]), int(line[8]), int(line[10]), line[11]))
        cdict[c2].append((line[12] + ":::end:::" + str(i + 1), int(line[3]), int(line[8]), int(line[10]), line[11]))
        if line[12] not in contig_coords:
            contig_coords[line[12]] = []
        contig_coords[line[12]].append(int(line[2]))
        contig_coords[line[12]].append(int(line[3]))
    coords[ref] = cdict




ref_segments = dict()
ref_segments_dict = {}
for ref, cs in coords.items():
    ref_segments_dict[ref] = []
for ref, cs in coords.items():
    scs = sorted(cs.items(), key=lambda z: z[0])
    for i in range(1, len(scs)):
        refname = scs[i][1][0][4]
        c1 = int(scs[i - 1][0])
        c2 = int(scs[i][0])
        e1 = (refname, c1)
        e2 = (refname, c2)
        if not e1 in ref_segments:
            ref_segments[e1] = []
        if not e2 in ref_segments:
            ref_segments[e2] = []
        ref_segments[e1].append((e1, e2))
        ref_segments[e2].append((e1, e2))
        ref_segments_dict[ref].append((e1, e2))


print ref_segments_dict, "REF_SEGMENTS_DICT"



for ref, cs in coords.items():
    print ref
    for u in sorted(cs.items(), key=lambda z: z[0]):
        print u
    print




working_contigs = set()
for ref, cs in coords.items():
    for u in cs.items():
        for v in u[1]:
            working_contigs.add(v[0].split(":::")[0])


# I need to determine which reference piece appears in which contigs
# sometimes, it may be misleading. We have to see if the piece is between ones piece start and end
# then we include it into the corresponding contig

contigs_refpieces = {}
for contig in working_contigs:
    contigs_refpieces[contig] = []
for ref, cs in coords.items():
    sorted_by_global_ref = sorted(cs.items(), key=lambda z: z[0])
    bag = dict()
    for ppp in range(len(sorted_by_global_ref[0][1])):
        cc, st, on = sorted_by_global_ref[0][1][ppp][0].split(":::")
        bag[on] = (cc, st)
    for uu in range(1, len(sorted_by_global_ref)):
        processed = []
        u = sorted_by_global_ref[uu]
        for v in u[1]: # iterate through contigs
            contig, status, order = v[0].split(":::")
            for ooo, info in bag.items():
                cntg, stts = info
                if stts == "start":
                    contigs_refpieces[cntg].append(((ref, sorted_by_global_ref[uu - 1][0]), (ref, sorted_by_global_ref[uu][0])))
            processed.append((order, contig, status))
        print processed, "PROCESSED"
        for rdr, cntg, stts in processed:
            if stts == "start":
                bag[rdr] = (cntg, stts)
            elif stts == "end":
                del bag[rdr]
        print bag, "BAG"

for ppp, pppp in contigs_refpieces.items():
    print ppp, pppp, "THIS IS CONTIGS"




breakpoints = {}
for contig in working_contigs:
    breakpoints[contig] = {}


for ref, cs in coords.items():
    bag = dict()
    for u in sorted(cs.items(), key=lambda z: z[0]):
        global_ref_coord = int(u[0])
        for v in u[1]: # iterate through contigs
            contig, status, order = v[0].split(":::")
            local_coord = int(v[1])
            contig_len = int(v[2])
            orient = int(v[3])
            refname = v[4]
            if status == "start":
                bag[order] = (global_ref_coord, contig, status, local_coord, orient, refname)
            for o, info in bag.items():
                grc, c, s, l, oo, r = info
                cur_point_local_delta = oo * (global_ref_coord - grc) # should be positive
                cur_point_local = l + cur_point_local_delta
                breakpoints[c][(refname, global_ref_coord)] = max(cur_point_local, 1)
            if status == "end":
                del bag[order]




segments_artificial = {}



for u, v in breakpoints.items():
    cdict = {}

    for k, t in v.items():
        if t not in cdict:
            cdict[t] = []
        cdict[t].append(k)


    good_segments = contigs_refpieces[u]

    print u, sorted(cdict.items(), key=lambda z: z[0]), len(good_segments), good_segments


    # let's see which points remain after we filter out small pieces (take for now 200 bp)


    sorted_points_on_contig = sorted(cdict.items(), key=lambda z: z[0])
    remaining_points = set()
    remaining_segments = set()
    for g in range(1, len(sorted_points_on_contig)):
        if sorted_points_on_contig[g][0] - sorted_points_on_contig[g - 1][0] > smallest_contig_length:
            remaining_points.add(g - 1)
            remaining_points.add(g)
            remaining_segments.add((g - 1, g))

    remaining_points = sorted(list(remaining_points))
    print remaining_points, "REMAINING POINTS"
    print remaining_segments, "REMAINING SEGMENTS"
    # now for each segment we have to find out it small covering sub-segments from the contig

    sorted_points_on_contig = sorted(cdict.items(), key=lambda z: z[0])
    for p1, p2 in good_segments:
        print "SEARCHING FOR", p1, p2
        gaps = []
        i1 = None
        i2 = None
        for i in range(len(sorted_points_on_contig)):
            lcoord, pts = sorted_points_on_contig[i]
            if p1 in pts:
                i1 = i
                break
        for i in range(len(sorted_points_on_contig)):
            lcoord, pts = sorted_points_on_contig[i]
            if p2 in pts:
                i2 = i
                break
        si1, si2 = tuple(sorted([i1, i2]))
        print si1, si2, "SI1 and SI2"

        valid_segments = set()
        for mmm in range(si1 + 1, si2 + 1):
            if (mmm - 1, mmm) in remaining_segments:
                valid_segments.add((mmm - 1, mmm))
        print valid_segments, "VS"
        sorted_valid_segments = sorted(list(valid_segments), key = lambda zzz: zzz[0])
        if not sorted_valid_segments:
            continue
        print sorted_valid_segments, "SVS"
        good_valid_segments = [[sorted_valid_segments[0]]]
        for mmm in range(1, len(sorted_valid_segments)):
            if sorted_valid_segments[mmm][0] == good_valid_segments[-1][-1][-1]:
                good_valid_segments[-1].append(sorted_valid_segments[mmm])
            else:
                good_valid_segments.append([sorted_valid_segments[mmm]])
        print good_valid_segments, "GVS"


        allsegm = set()
        for i in range(si1, si2):
            allsegm.add((i, i + 1))
        gapsegm = sorted(list(allsegm - set(sorted_valid_segments)), key=lambda zzz: zzz[0])
        print gapsegm, "GAPSEGM"

        if gapsegm:
            gap_segments = [[gapsegm[0]]]
            for mmm in range(1, len(gapsegm)):
                if gapsegm[mmm][0] == gap_segments[-1][-1][-1]:
                    gap_segments[-1].append(gapsegm[mmm])
                else:
                    gap_segments.append([gapsegm[mmm]])
            gap_segments = map(lambda zzz: (zzz[0][0], zzz[-1][1]), gap_segments)
            gaps2 = map(lambda zzz: sorted_points_on_contig[zzz[1]][0] - sorted_points_on_contig[zzz[0]][0], gap_segments)
            print gap_segments, gaps2, "GAP SEGMENTS"
        else:
            gap_segments = []
            gaps2 = []
            print gap_segments, "GAP SEGMENTS"


        gapcontig_dict = {}
        for xxx in gapsegm:
            gapcontig_dict[xxx] = "GAP"
        for xxx in sorted_valid_segments:
            gapcontig_dict[xxx] = "CONTIG"
        gapcontig_list = []
        for xxx in sorted(gapcontig_dict.keys(), key=lambda zzz: zzz[0]):
            if gapcontig_list == []:
                gapcontig_list.append(gapcontig_dict[xxx])
            elif gapcontig_dict[xxx] == gapcontig_list[-1]:
                pass
            else:
                gapcontig_list.append(gapcontig_dict[xxx])
        print gapcontig_list, "GAPCONTIG_LIST"

        good_points = []
        for xxx in good_valid_segments:
            for xxxx in xxx:
                good_points.append(xxxx)
        #good_points = map(lambda mmm: (mmm[0][0], mmm[-1][1]), good_valid_segments)
        #print good_points, "GPS"


        final_subsegments = []
        for gp1, gp2 in good_points:
            it1 = sorted_points_on_contig[gp1]
            it2 = sorted_points_on_contig[gp2]
            print "GP GP", gp1, gp2, it1, it2, p1, p2
            if i1 < i2:
                final_subsegments.append("%s:%s-%s" % (u, it1[0], it2[0]))
            else:
                final_subsegments.append("%s:%s-%s" % (u, it2[0], it1[0]))
        print "%s\t%s\t%s-%s" % ("IGOR_MANDRIC", final_subsegments, p1, p2)


        gaps22 = iter(gaps2)
        final_subsegments2 = iter(final_subsegments)
        final_chain = []
        for matter in gapcontig_list:
            if matter == "GAP":
                final_chain.append(gaps22.next())
            elif matter == "CONTIG":
                final_chain.append(final_subsegments2.next())
        if i1 > i2:
            final_chain = list(reversed(final_chain))
        print "FINAL CHAIN IS THIS", final_chain


        if (p1, p2) not in segments_artificial:
            segments_artificial[(p1, p2)] = {}
        segments_artificial[(p1, p2)][u] = final_chain



print "-----------------------------------------------------------------"



scafdict = {}
for aaaa, bbbb in ref_segments_dict.items():
    scafdict[aaaa] = []


for aaaa, bbbb in ref_segments_dict.items():
    print aaaa, "THIS IS REFERENCE"
    for xxx in bbbb:
        print xxx, "\t", segments_artificial.get(xxx)
        piece = segments_artificial.get(xxx)
        if piece is None:
            scafdict[aaaa].append(xxx[1][1] - xxx[0][1])
        else:
            alphabetically = sorted(piece.keys())[0]
            scafdict[aaaa].extend(piece[alphabetically])


for aaaa, bbbb in scafdict.items():
    print "REFERENCE", aaaa
    for xxx in bbbb:
        print xxx


# now combine the gaps

scafdict2 = {}

for aaaa, bbbb in scafdict.items():
    i = 0
    while i < len(bbbb) and type(bbbb[i]) == int:
        i += 1
        pass
    if i >= len(bbbb):
        continue
    scaffold = [bbbb[i]]
    for xxx in bbbb[i + 1:]:
        if type(xxx) == int and type(scaffold[-1]) == int:
            scaffold[-1] += xxx
        elif type(xxx) == str and type(scaffold[-1]) == str:
            scaffold.append(0)
            scaffold.append(xxx)
        else:
            scaffold.append(xxx)
    if type(scaffold[-1]) == int:
        scaffold[-1] = 0
    else:
        scaffold.append(0)
    scafdict2[aaaa] = scaffold


for aaaa, bbbb in scafdict2.items():
    print "REFERENCE MMMMM", aaaa
    for xxx in bbbb:
        print xxx


final_scaffolds = []
fasta = {}
used_or_not = dict()
counter = 1
for aaaa, bbbb in scafdict2.items():
    print len(bbbb), "AAAAAAAAAAAAAAAAAAAAAAAA"
    scaffold = []
    for i in range(0, len(bbbb), 2):
        contigname, coords = bbbb[i].split(":")
        coord1, coord2 = coords.split("-")
        cc1, cc2 = tuple(sorted([int(coord1), int(coord2)]))
        toused = ":".join(map(str, [contigname, cc1, cc2]))
        sequence = sw[contigname][cc1 - 1: cc2]
        if toused not in used_or_not:
            name = "contig_%s" % counter
            used_or_not[toused] = name
            counter += 1
        else:
            name = used_or_not[toused]
        print name, bbbb[i]
        fasta[name] = sequence
        if int(coord1) < int(coord2):
            orientation = 1
        else:
            orientation = -1
        scaffold.append("%s:::%s:::%s:::%s\n" % (name, orientation, len(sequence), bbbb[i + 1]))
        print len(scaffold), "JJJKJK"
    final_scaffolds.append(scaffold)


with open(prefix + ".scaf", "w") as f:
    for i, scaffold in enumerate(final_scaffolds):
        f.write(">scaffold_%s\n" % (i + 1))
        for contig in scaffold:
            f.write(contig)


with open(prefix + ".fasta", "w") as f:
    for xxx in range(len(fasta)):
        f.write(">contig_%s\n" % (xxx + 1))
        f.write("%s\n" % fasta["contig_%s" % (xxx + 1)])



















