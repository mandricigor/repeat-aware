


import sys
import logging
#import numpy as np
import time
import cplex
import math
#import networkx as nx
from cplex.exceptions import CplexSolverError
from itertools import tee, izip
import argparse

def op(ori):
    if ori == "+":
        return "-"
    else:
        return "+"

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)




# toy example
#scaffold = [['M', 'A', 'B'], ['C', 'D', 'G'], ['E', 'H', 'I'], ['A', 'B', 'C', 'D', 'E', 'F'], ['G', 'F', 'J', 'K', 'A']]
#scaf_orient = [["+", "+", "+"], ["+", "+", "+"], ["+", "+", "+"], ["+", "+", "+", "+", "+", "+"], ["+", "+", "+", "+", "-"]]

#reference = [['M', 'A', 'B', 'C', 'D', 'G', 'E', 'H', 'I', 'A', 'B', 'C', 'G', 'F', 'J', 'K', 'D', 'E', 'F', 'Z']]
#ref_orient = [["+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "-"]]



def parse_scaf_file(scaffile):
    with open(scaffile) as f:
        a = f.readlines()
    a = "".join(a)
    a = a.split(">")[1:]
    a = map(lambda x: x.split("\n")[1:-1], a)
    answer = []
    orient = []
    distances = []
    lengths_main = []
    for scaf in a:
        scaf2 = []
        ori = []
        dist = []
        lengths = []
        for contig in scaf:
            contig = contig.split(":::")
            name = contig[0].split("_")[1]
            gap = int(contig[3])
            length = int(contig[2])
            dist.append(gap)
            lengths.append(length)
            scaf2.append(name)
            if contig[1] == "1":
                ori.append("+")
            else:
                ori.append("-")
        answer.append(scaf2)
        orient.append(ori)
        del dist[-1]
        distances.append(dist)
        lengths_main.append(lengths)
    return answer, orient, distances, lengths_main





parser = argparse.ArgumentParser(description='Repeat-aware evaluation framework of Igor Mandric')
parser.add_argument('--outscaf', dest="outscaf", type=str, required=True, help="output scaffolding")
parser.add_argument('--refscaf', dest="refscaf", type=str, required=True, help="reference scaffolding")
parser.add_argument('--dist', dest="distance_threshold", type=int, required=False, default=1000, help="links for which estimated gap differs from the true value by more than DISTANCE_THRESHOLD bp are considered wrong (default 1000 bp)")
parser.add_argument('--noskip', dest="noskip", type=int, required=False, help="if the argument is present than do not consider links skipping contigs as wrong if the contigs are at most NOSKIP bp apart on the reference, otherwise links skipping contigs are wrong")


args = parser.parse_args()

scaffold, scaf_orient, sdistances, slengths = parse_scaf_file(args.outscaf)
reference, ref_orient, rdistances, rlengths = parse_scaf_file(args.refscaf)
distance_threshold = args.distance_threshold
if args.noskip is None:
    noskip = False
    noskip_distance = 0
else:
    noskip = True
    noskip_distance = args.noskip



# paper example
#scaffold = [["f", "a", "b", "c", "d", "e", "d"], ["g"]]
#scaf_orient = [["+", "+", "+", "-", "+", "+", "+"], ["-"]]

#reference = [["a", "b", "c", "b", "d", "e"], ["f", "g"]]
#ref_orient = [["+", "+", "+", "-", "+", "+"], ["+", "+"]]




# dictionary containing relationship: contig -> [reference]
refdict = {}
for count, ref in enumerate(reference):
    for r in ref:
        if r not in refdict:
            refdict[r] = set()
        refdict[r].add(count)








allscaf = [item for sublist in scaffold for item in sublist]
allref = [item for sublist in reference for item in sublist]


sorient = [item for sublist in scaf_orient for item in sublist]
rorient = [item for sublist in ref_orient for item in sublist]





cons1 = []
cons2 = []
condict = {}


ref_distances = {}
i = 0
for k, refer in enumerate(reference):
    for j in range(len(refer) - 1):
        ref_distances[(i, i + 1)] = rdistances[k][j]
        i += 1
    i += 1


scaf_distances = {}
i = 0
for k, sc in enumerate(scaffold):
    for j in range(len(sc) - 1):
        scaf_distances[(i, i + 1)] = sdistances[k][j]
        i += 1
    i += 1




if noskip == True:

    r_cont_number = sum(map(len, reference))

    i = 0
    for k in range(r_cont_number - 1):
        #print i, i + 1, "MKKK"
        if (i, i + 1) not in ref_distances:
            ref_distances[(i, i + 1)] = 100000000
        i += 1


    allrlens = []
    for r in rlengths:
        allrlens.extend(r)


    ref_all_distances = {}
    # computed distances between contigs on the reference
    for k in range(r_cont_number - 1):
        overall_distance = ref_distances[(k, k + 1)] + allrlens[k + 1]
        #print overall_distance, k, "KJKJKJK"
        for m in range(k + 2, r_cont_number):
            overall_distance += ref_distances[(m - 1, m)]
            #print k, m, overall_distance, "JJJJJ"
            ref_all_distances[(k, m)] = overall_distance
            overall_distance += allrlens[m]


    print ref_all_distances, "REF_ALL_DISTANCES"

    for xxx in ref_distances:
        ref_all_distances[xxx] = ref_distances[xxx]

else:
    ref_all_distances = ref_distances









scaf_distances = {}


variables = []
i = 0
for ref in reference:
    for r in ref:
        j = 0
        scafn = 0
        c1 = []
        for scaf in scaffold:
            k = 0
            for sc in scaf:
                if r == sc:
                    pos = "%s_%s" % (scafn, j)
                    variables.append((i, pos))
                    c1.append((i, pos))
                    p = condict.get(pos, [])
                    p.append(i)
                    condict[pos] = p
                # add distances here to the dictionary
                if k + 1 < len(scaffold[scafn]):
                    scaf_distances[("%s_%s" % (scafn, j), "%s_%s" % (scafn, j + 1))] = sdistances[scafn][k]
                k += 1
                j += 1
            scafn += 1
        i += 1
        cons1.append(c1)

#print cons1, "CONS1"



for x, y in condict.items():
    if len(y) > 1:
        cons = []
        for yy in y:
            cons.append((yy, x))
        cons2.append(cons)




# inter contig constraints
scaf_outer = set()
for scaf in scaffold:
    for a1, a2 in pairwise(scaf):
        scaf_outer.add(tuple(sorted([a1, a2])))


#print scaf_outer, "SCAF OUTER"



counts = {}
for scaf in scaffold:
    for sc in scaf:
        counts[sc] = 0


outer_indices = {}
outer_orient = {}
indices = []

for scaf in scaffold:
    for a in scaf:
        indices.append(counts[a])
        counts[a] += 1


j = 0
for scaf, ori in zip(scaffold, scaf_orient):
    for i in range(len(scaf) - 1):
        a, b = scaf[i], scaf[i + 1]
        o1, o2 = ori[i], ori[i + 1]
        link = tuple(sorted([a, b]))
        if link not in outer_indices:
            outer_indices[link] = []
        if link not in outer_orient:
            outer_orient[link] = []
        if link == (a, b):
            outer_indices[link].append([indices[j], indices[j + 1]])
            outer_orient[link].append([o1, o2])
        else:
            outer_indices[link].append([indices[j + 1], indices[j]])
            outer_orient[link].append([op(o2), op(o1)])
        j += 1
    j += 1





cons_pairs_outer = {}
cons_pairs_orient = {}

"""
j = 0
for ref, ori in zip(reference, ref_orient):
    for i in range(len(ref) - 1):
        link = tuple(sorted([ref[i], ref[i + 1]]))
        o1, o2 = ori[i], ori[i + 1]
        if link in scaf_outer:
            if link not in cons_pairs_outer:
                cons_pairs_outer[link] = []
            if link not in cons_pairs_orient:
                cons_pairs_orient[link] = []
            if link == (ref[i], ref[i + 1]):
                cons_pairs_outer[link].append([cons1[j], cons1[j + 1]])
                cons_pairs_orient[link].append([o1, o2])
            else:
                cons_pairs_outer[link].append([cons1[j + 1], cons1[j]])
                cons_pairs_orient[link].append([op(o2), op(o1)])
        j += 1
    j += 1
"""

for crefs, di in ref_all_distances.items():
    if (not noskip) or (noskip and di < noskip_distance):
        cref1, cref2 = crefs
        link = tuple(sorted([allref[cref1], allref[cref2]]))
        #print link, "LINKLLLLLLLLLLLLLLLLLLL"
        if link not in outer_indices:
            continue
        o1, o2 = rorient[cref1], rorient[cref2]
        if link not in cons_pairs_outer:
            cons_pairs_outer[link] = []
        if link not in cons_pairs_orient:
            cons_pairs_orient[link] = []
        if link == (allref[cref1], allref[cref2]):
            cons_pairs_outer[link].append([cons1[cref1], cons1[cref2]])
            cons_pairs_orient[link].append([o1, o2])
        else:
            cons_pairs_outer[link].append([cons1[cref2], cons1[cref1]])
            cons_pairs_orient[link].append([op(o2), op(o1)])
        






# now build the cplex problem

cpx = cplex.Cplex()
cpx.set_results_stream("solution.txt")


# add variables

for var in variables:
    Xij = "X#%s#%s" % var
    cpx.variables.add(lb = [0], ub = [1], types = ["B"], names = [Xij])


# constraints 1
for con in cons1:
    inds = map(lambda x: "X#%s#%s" % x, con)
    vals = [1] * len(inds)
    c = cplex.SparsePair(ind = inds, val = vals)
    if len(vals) == 1:
        sens = "L"
    else:
        sens = "L"
    cpx.linear_constraints.add( \
        lin_expr = [c],\
        senses = [sens],\
        rhs = [1],\
        names = ['pair']\
    )



# constraints 2
for con in cons2:
    inds = map(lambda x: "X#%s#%s" % x, con)
    vals = [1] * len(inds)
    c = cplex.SparsePair(ind = inds, val = vals)
    cpx.linear_constraints.add( \
        lin_expr = [c],\
        senses = ["L"],\
        rhs = [1],\
        names = ['pair']\
    )

#print ref_distances
#print scaf_distances



#print cons_pairs_outer, "HHHHHHHHHHHHHHHHHHHHHHHHHHHHH"



allzeds = []
zedweights = []

wrong_distance_links = set()
xz_corr = {} # correspondence between x variables and z
used_x = set()
z = 0
for link in cons_pairs_outer:
    #print link, "THIS IS LINK"
    conses = cons_pairs_outer[link]
    conses_orient = cons_pairs_orient[link]
    cind = outer_indices[link]
    cind_orient = outer_orient[link]
    #zeds = []
    for cons, cori in zip(conses, conses_orient):
        #print cons, cori, "CONS CORI"
        zz = []
        ww = []
        for index, sori in zip(cind, cind_orient):
            rgap = max(0, ref_all_distances[tuple(sorted([cons[0][index[0]][0], cons[1][index[1]][0]]))])
            sgap = max(0, scaf_distances[tuple(sorted([cons[0][index[0]][1], cons[1][index[1]][1]], key=lambda mm: int(mm.split("_")[1])))])
            #print sgap, rgap, cons, cind, index
            if abs(rgap - sgap) > distance_threshold:
                #print "ADDING WRONG", tuple(sorted([int(cons[0][index[0]][0]), int(cons[1][index[1]][0])]))
                #wrong_distance_links.add(tuple(sorted([int(cons[0][index[0]][1].split("_")[1]), int(cons[1][index[1]][1].split("_")[1])])))
                wrong_distance_links.add(tuple(sorted([int(cons[0][index[0]][0]), int(cons[1][index[1]][0])])))
            elif cori == sori and abs(rgap - sgap) <= distance_threshold:
                #print "ADDING CORRECT"
                w = 1
                ww.append(w)
                cpx.variables.add(lb = [0], ub = [1], types = ["B"], names = ["Z#%s" % z])
                inds = ["X#%s#%s" % cons[0][index[0]], "X#%s#%s" % cons[1][index[1]], "Z#%s" % z]
                #print inds
                used_x.add("X#%s#%s" % cons[0][index[0]])
                used_x.add("X#%s#%s" % cons[1][index[1]])
                xz_corr["Z#%s" % z] = ("X#%s#%s" % cons[0][index[0]], "X#%s#%s" % cons[1][index[1]])
                vals = [1, 1, -2]
                c = cplex.SparsePair(ind = inds, val = vals)
                cpx.linear_constraints.add( \
                    lin_expr = [c],\
                    senses = ["G"],\
                    rhs = [0],\
                    names = ['pair']\
                )
                zz.append(["Z#%s" % z])
                z += 1
            #else:
                #print cori, sori, "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
        for zzz in zz:
            allzeds.extend(zzz)    
        zedweights.extend(ww)
        #print zz, "THIS  IS ZZ"
        if len(zz) > 1:
            for zzz in zip(*zz):
                #print zzz, "JJJJJJJJJJJJJJJJJJJJJJJJJJ"
                c = cplex.SparsePair(ind = zzz, val = [1] * len(zz))
                cpx.linear_constraints.add( \
                    lin_expr = [c],\
                    senses = ["L"],\
                    rhs = [1],\
                    names = ['pair']\
                )

"""
        zeds.append(zz)
        zedweights.extend(ww)
    for zzz in zeds:
        allzeds.extend(zzz)
    print zeds, "ZEDS"
    if len(zeds) > 1:
        for zzz in zip(*zeds):
            c = cplex.SparsePair(ind = zzz, val = [1] * len(zeds))
            cpx.linear_constraints.add( \
                lin_expr = [c],\
                senses = ["L"],\
                rhs = [1],\
                names = ['pair']\
            )
"""

#print allzeds, "ALLZEDS"


for zed, w in zip(allzeds, zedweights):
    cpx.objective.set_linear(zed, w)



cpx.objective.set_sense(cpx.objective.sense.maximize)
cpx.set_problem_type(cpx.problem_type.MILP)
cpx.write("program.txt", filetype="lp")
start_time = cpx.get_time()
cpx.solve()
elapsed = cpx.get_time() - start_time

answer_dict = {}

for x in cpx.variables.get_names():
    if "X" in x:
        x = x.split("#")[1]
        answer_dict[x] = []


# --------------------------- processing ------------------

#for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
#    print x, y, "JJJJJJJJJJJJJ"

scafs = []
for scaf in scaffold:
    for sc in scaf:
        scafs.append("NONE")

for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
    if y > 0.5 and "X" in x:
        x = x.split("#")
        a = int(x[1])
        b = int(x[2].split("_")[1])
        scafs[b] = a

good_links = []
for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
    #if "Z" in x:
    #    print x, y, "MMM"
    if y > 0.5 and "Z" in x:
        # collect correct links
        fl, sl = xz_corr[x]
        fl = fl.split("#")
        b1 = int(fl[2].split("_")[1])
        sl = sl.split("#")
        b2 = int(sl[2].split("_")[1])
        #print b1, b2, "GOOD LINK"
        correct = tuple(sorted([b1, b2]))
        good_links.append(correct)



j = 0
scaf_answer = []
for scaf in scaffold:
    sa = []
    for sc in scaf:
        sa.append(scafs[j])
        j += 1
    scaf_answer.append(sa)

# now fill out those which were not assigned to anything
assigned = {}
for i, sc in enumerate(scaf_answer):
    for j, s in enumerate(sc):
        if s != "NONE":
            assigned[s] = (i, j)


newly_assigned = {}
for i, x, y in zip(range(len(scaf_answer)), scaf_answer, scaffold):
    for j, xx, yy in zip(range(len(x)), x, y):
        if xx == "NONE":
            indices = []
            for k in range(len(allref)):
                if allref[k] == yy:
                    indices.append(k)
            for ind in indices:
                if ind in assigned:
                    continue
                else:
                    assigned[ind] = (i, j)
                    newly_assigned[ind] = (i, j)
                    break



for x, y in newly_assigned.items():
    u, v = y
    scaf_answer[u][v] = x


#print scaf_answer

# after assignment, if there is still something 
scaf_answer_human = []
asc = iter(allscaf)
for scaf in scaf_answer:
    myscaf = []
    for sc in scaf:
        myscaf.append(asc.next())
    scaf_answer_human.append(myscaf)

#print scaf_answer_human
#print len(scaf_answer_human)
#print map(len, scaf_answer_human)

#print good_links
wrong_links = []
corrected_chains = []
position = 0
for scaffold in scaf_answer:
    correct_subscaffold = [[scaffold[0]]]
    for j in range(1, len(scaffold)):
        if (position, position + 1) in good_links:
            print "CORRECT LINK:", scaffold[j - 1], scaffold[j], allscaf[position], allscaf[position + 1]
            correct_subscaffold[-1].append(scaffold[j])
        else:
            correct_subscaffold.append([scaffold[j]])
            wrong_links.append((position, position + 1, scaffold[j - 1], scaffold[j]))
        position += 1
    position += 1
    corrected_chains.extend(correct_subscaffold)

#print corrected_chains, "CORRECTED CHAINS"

# compute correct links
#correct_links = []
#for chain in corrected_chains:
#    if len(chain) > 1:
#        for i in range(1, len(chain)):
#            correct_links.append((allref[chain[i - 1]], allref[chain[i]]))
        
#for u, v in correct_links:
#    print u, v

#print wrong_links, "WRONG LINKS"

#print corrected_chains


scaf_answer_human = []
scaf_answer_nonhuman = []
asc = iter(allscaf)
asc2 = iter(range(len(allscaf)))
for scaf in corrected_chains:
    myscaf = []
    myscaf2 = []
    for sc in scaf:
        myscaf.append(asc.next())
        myscaf2.append(asc2.next())
    scaf_answer_human.append(myscaf)
    scaf_answer_nonhuman.append(myscaf2)


scaf_answer_nonhuman2 = []
for uuu, vvv in zip(scaf_answer_nonhuman, corrected_chains):
    if vvv != ['NONE']:
        scaf_answer_nonhuman2.append(uuu)

corrected_chains = [x for x in corrected_chains if x != ['NONE']]

#print corrected_chains
#print scaf_answer_nonhuman2, "NONHUMAN"

# ---------------- code for computing corrected n50

with open(args.outscaf) as f:
    lines = f.readlines()
lines = map(lambda x: x.strip(), lines)
lines = [x for x in lines if ">" not in x]

corlines = []
for cchain, vvv in zip(scaf_answer_nonhuman2, corrected_chains):
    corchain = []
    for number in cchain:
        corchain.append(lines[number])
    #if vvv[0] > vvv[-1]:
    #    corlines.append(list(reversed(corchain)))
    #else:
    corlines.append(corchain)


corlens = []
for cline in corlines:
    l1 = map(lambda x: int(x.split(":::")[2]), cline)
    l2 = map(lambda x: int(x.split(":::")[3]), cline)
    corlens.append(sum(l1) + sum(l2[:-1]))


sumalllens = sum(corlens) * 1.0 / 2
scorlens = sorted(corlens, reverse=True)
cursum = 0
i = 0
while (cursum < sumalllens):
    cursum += scorlens[i]
    i += 1

#print "CORRECTED N50:", scorlens[i - 1]
corn50 = scorlens[i - 1]




# -------------------------------------------------


#print scaf_answer_human

# classify wrong links
wrong_copy = 0
jumping = 0
wrong_ord_ori = 0
wrong_ref = 0
jumping_wrong_ord_ori = 0
wrong_distance = 0


#print wrong_distance_links

for aa, bb, cc, dd in wrong_links:

    print aa, bb, cc, dd, allscaf[aa], allscaf[bb] # enumerate all wrong links
    is_jumping = False
    is_wrong_ord_ori = False
    is_wrong_copy = False
    is_wrong_ref = False


    if cc != "NONE" and dd != "NONE":
        o1, o2, o3, o4 = sorient[int(aa)], sorient[int(bb)], rorient[int(cc)], rorient[int(dd)]

        is_wrong_ref = (refdict[allref[int(cc)]] & refdict[allref[int(dd)]] == set())
        if is_wrong_ref:
            print aa, bb, cc, dd, "WRONG_REF", refdict[allref[int(cc)]], refdict[allref[int(dd)]], refdict[allref[int(cc)]] & refdict[allref[int(dd)]], refdict[allref[int(cc)]] & refdict[allref[int(dd)]] == set(), is_wrong_ref, allscaf[aa], allscaf[bb]
            wrong_ref += 1
            #print aa, bb, cc, dd, "MMMMMMMMMMM", allscaf[aa], allscaf[bb]
            continue
        is_wrong_copy = False
        is_jumping = abs(int(cc) - int(dd)) > 1
        #print is_jumping, allscaf[aa], allscaf[bb], "IS JUMPING"
        #print aa < bb, cc < dd, (o1, o2 == o3, o4), (o1, o2 == op(o4), op(o3)), "KKK"
        if not (((aa < bb) and (cc < dd) and ((o1, o2) == (o3, o4))) or ((((aa < bb) and (cc > dd)) or ((aa > bb) and (cc < dd))) and ((o1, o2) == (op(o3), op(o4))))):
            is_wrong_ord_ori = True
            print aa, bb, cc, dd, "WRONG ORDORI", allscaf[aa], allscaf[bb], o1, o2, o3, o4


        #print aa, bb, cc, dd, is_wrong_ord_ori, tuple(sorted([cc, dd])) in wrong_distance_links, o1, o2, o3, o4
        if tuple(sorted([cc, dd])) in wrong_distance_links and is_wrong_ord_ori == False:
            wrong_distance += 1
            continue

    else:
        print aa, bb, cc, dd, "MMMMMMMMMMMMMMM", allscaf[aa], allscaf[bb], "WRONG COPY"
        is_jumping = False
        is_wrong_ord_ori = False
        is_wrong_ref = False
        is_wrong_copy = True

    #print is_wrong_copy, is_wrong_ref, is_jumping, is_wrong_ord_ori, wrong_distance

    if is_wrong_copy:
        wrong_copy += 1
    #elif is_wrong_ref:
    #    wrong_ref += 1
    elif is_jumping and is_wrong_ord_ori:
        jumping_wrong_ord_ori += 1
    elif is_jumping:
        jumping += 1
    elif is_wrong_ord_ori:
        wrong_ord_ori += 1

report = []

try:
    tool = sys.argv[1].split("/")[1]
except Exception:
    tool = "TOOL"
report.append(tool)

# put info about species, sim, len
try:
    species, sim, length, suf = sys.argv[2].split("/")[-1].split(".")
except Exception:
    species, sim, length, suf = "A", "B", "C", "D"
report.extend([species, sim, length])

#print "COPY:" , wrong_copy
#print "REF:", wrong_ref
#print "GAP:", wrong_distance
#print "JUMPING:", jumping
#print "ORDER_ORIENTATION:", wrong_ord_ori
#print "J+ORD_ORI", jumping_wrong_ord_ori

report.extend([wrong_copy, wrong_ref, wrong_distance, jumping, wrong_ord_ori, jumping_wrong_ord_ori])

corlens = map(len, corrected_chains)
#print "MEAN CHAIN LENGTH:", sum(corlens) * 1.0 / len(corlens)
#print "MAX CHAIN LENGTH:", max(corlens)

report.extend([float("%.2f" % (sum(corlens) * 1.0 / len(corlens))), max(corlens)])

# compute n50 (contig chains)
sumalllens = sum(corlens) * 1.0 / 2
scorlens = sorted(corlens, reverse=True)
cursum = 0
i = 0
while (cursum < sumalllens):
    cursum += scorlens[i]
    i += 1

#print "CORRECTED CHAIN N50:", scorlens[i - 1]
report.append(scorlens[i - 1])
report.append(corn50)

correct_links = int(cpx.solution.get_objective_value())

#print "CORRECT LINKS:", correct_links
report.append(correct_links)



inferred_links =  sum(map(lambda x: x - 1, map(len, scaf_answer)))
report.append(inferred_links)

ref_links = sum(map(lambda x: x - 1, map(len, reference)))


sens = correct_links * 1.0 / ref_links

if inferred_links > 0:
    ppv = correct_links * 1.0 / inferred_links
else:
    ppv = 100.0
#print "%.2f" % sens, "%.2f" % ppv

report.extend([float("%.2f" % sens), float("%.2f" % ppv)])


print "TOTAL LINKS", ref_links
print "CORRECT LINKS", correct_links
print "INFERRED LINKS", inferred_links


sens = float("%.2f" % sens)
ppv = float("%.2f" % ppv)

if sens == 0 or ppv == 0:
    fscore = 0
else:
    fscore = 2 * sens * ppv / (sens + ppv)
report.append(float("%.2f" % fscore))

report.append("%.2f" % elapsed)
print ",".join(map(str, report))


"""
COPY: 0
REF: 1
GAP: 0
JUMPING: 0
ORDER_ORIENTATION: 0
J+ORD_ORI 47
MEAN CHAIN LENGTH: 1.18367346939
MAX CHAIN LENGTH: 7
CHAIN N50: 1
CORRECT LINKS: 108
0.13 0.69
"""

