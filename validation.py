


import sys
import logging
#import numpy as np
import cplex
import math
#import networkx as nx
from cplex.exceptions import CplexSolverError
from itertools import tee, izip


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
    for scaf in a:
        scaf2 = []
        ori = []
        dist = []
        for contig in scaf:
            contig = contig.split(":::")
            name = contig[0].split("_")[1]
            gap = int(contig[3])
            dist.append(gap)
            scaf2.append(name)
            if contig[1] == "1":
                ori.append("+")
            else:
                ori.append("-")
        answer.append(scaf2)
        orient.append(ori)
        del dist[-1]
        distances.append(dist)
    return answer, orient, distances







scaffold, scaf_orient, sdistances = parse_scaf_file(sys.argv[1])
reference, ref_orient, rdistances = parse_scaf_file(sys.argv[2])
if len(sys.argv) < 4:
    distance_threshold = 1000
else:
    distance_threshold = float(sys.argv[3])




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
#allref_distances = [item for sublist in rdistances for item in sublist]
i = 0
for k, refer in enumerate(reference):
    for j in range(len(refer) - 1):
        ref_distances[(i, i + 1)] = rdistances[k][j]
        i += 1
    i += 1

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



allzeds = []
zedweights = []

wrong_distance_links = set()
xz_corr = {} # correspondence between x variables and z
used_x = set()
z = 0
for link in cons_pairs_outer:
    conses = cons_pairs_outer[link]
    conses_orient = cons_pairs_orient[link]
    cind = outer_indices[link]
    cind_orient = outer_orient[link]
    zeds = []
    for cons, cori in zip(conses, conses_orient):
        zz = []
        ww = []
        for index, sori in zip(cind, cind_orient):
            rgap = max(0, ref_distances[tuple(sorted([cons[0][index[0]][0], cons[1][index[1]][0]]))])
            sgap = max(0, scaf_distances[tuple(sorted([cons[0][index[0]][1], cons[1][index[1]][1]], key=lambda mm: int(mm.split("_")[1])))])
            if abs(rgap - sgap) > distance_threshold:
                wrong_distance_links.add(tuple(sorted([int(cons[0][index[0]][1].split("_")[1]), int(cons[1][index[1]][1].split("_")[1])])))
            elif cori == sori and abs(rgap - sgap) <= distance_threshold:
                w = 1
                ww.append(w)
                cpx.variables.add(lb = [0], ub = [1], types = ["B"], names = ["Z#%s" % z])
                inds = ["X#%s#%s" % cons[0][index[0]], "X#%s#%s" % cons[1][index[1]], "Z#%s" % z]
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
                zz.append("Z#%s" % z)
                z += 1
        zeds.append(zz)
        zedweights.extend(ww)
    for zzz in zeds:
        allzeds.extend(zzz)
    if len(zeds) > 1:
        for zzz in zip(*zeds):
            c = cplex.SparsePair(ind = zzz, val = [1] * len(zeds))
            cpx.linear_constraints.add( \
                lin_expr = [c],\
                senses = ["L"],\
                rhs = [1],\
                names = ['pair']\
            )




for zed, w in zip(allzeds, zedweights):
    cpx.objective.set_linear(zed, w)



cpx.objective.set_sense(cpx.objective.sense.maximize)
cpx.set_problem_type(cpx.problem_type.MILP)
cpx.write("program.txt", filetype="lp")
cpx.solve()

answer_dict = {}
for x in cpx.variables.get_names():
    if "X" in x:
        x = x.split("#")[1]
        answer_dict[x] = []


# --------------------------- processing ------------------

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
    if y > 0.5 and "Z" in x:
        # collect correct links
        fl, sl = xz_corr[x]
        fl = fl.split("#")
        b1 = int(fl[2].split("_")[1])
        sl = sl.split("#")
        b2 = int(sl[2].split("_")[1])
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


# after assignment, if there is still something 
scaf_answer_human = []
asc = iter(allscaf)
for scaf in scaf_answer:
    myscaf = []
    for sc in scaf:
        myscaf.append(asc.next())
    scaf_answer_human.append(myscaf)




wrong_links = []
corrected_chains = []
position = 0
for scaffold in scaf_answer:
    correct_subscaffold = [[scaffold[0]]]
    for j in range(1, len(scaffold)):
        if (position, position + 1) in good_links:
            correct_subscaffold[-1].append(scaffold[j])
        else:
            correct_subscaffold.append([scaffold[j]])
            wrong_links.append((position, position + 1, scaffold[j - 1], scaffold[j]))
        position += 1
    position += 1
    corrected_chains.extend(correct_subscaffold)

corrected_chains = [x for x in corrected_chains if x != ['NONE']]


scaf_answer_human = []
asc = iter(allscaf)
for scaf in corrected_chains:
    myscaf = []
    for sc in scaf:
        myscaf.append(asc.next())
    scaf_answer_human.append(myscaf)




# classify wrong links
wrong_copy = 0
jumping = 0
wrong_ord_ori = 0
wrong_ref = 0
jumping_wrong_ord_ori = 0
wrong_distance = 0



for aa, bb, cc, dd in wrong_links:
    is_jumping = False
    is_wrong_ord_ori = False
    is_wrong_copy = False
    is_wrong_ref = False

    if (aa, bb) in wrong_distance_links:
        wrong_distance += 1
        continue
    if cc != "NONE" and dd != "NONE":
        o1, o2, o3, o4 = sorient[int(aa)], sorient[int(bb)], rorient[int(cc)], rorient[int(dd)]
        is_jumping = abs(int(cc) - int(dd)) > 1
        if (o1, o2 != o3, o4) or (o1, o2 != op(o4), op(o3)):
            is_wrong_ord_ori = True
        is_wrong_copy = False
        is_wrong_ref = refdict[allref[int(cc)]] & refdict[allref[int(dd)]] == set()
    else:
        is_jumping = False
        is_wrong_ord_ori = False
        is_wrong_ref = False
        is_wrong_copy = True

    if is_wrong_copy:
        wrong_copy += 1
    elif is_wrong_ref:
        wrong_ref += 1
    elif is_jumping and is_wrong_ord_ori:
        jumping_wrong_ord_ori += 1
    elif is_jumping:
        jumping += 1
    elif is_wrong_ord_ori:
        wrong_ord_ori += 1

print "COPY:" , wrong_copy
print "REF:", wrong_ref
print "GAP:", wrong_distance
print "JUMPING:", jumping
print "ORDER_ORIENTATION:", wrong_ord_ori
print "J+ORD_ORI", jumping_wrong_ord_ori


corlens = map(len, corrected_chains)
print "MEAN CHAIN LENGTH:", sum(corlens) * 1.0 / len(corlens)
print "MAX CHAIN LENGTH:", max(corlens)


# compute n50 (contig chains)
sumalllens = sum(corlens) * 1.0 / 2
scorlens = sorted(corlens, reverse=True)
cursum = 0
i = 0
while (cursum < sumalllens):
    cursum += scorlens[i]
    i += 1

print "CHAIN N50:", scorlens[i - 1]


correct_links = int(cpx.solution.get_objective_value())

print "CORRECT LINKS:", correct_links

inferred_links =  sum(map(lambda x: x - 1, map(len, scaf_answer)))

ref_links = sum(map(lambda x: x - 1, map(len, reference)))


sens = correct_links * 1.0 / ref_links

if inferred_links > 0:
    ppv = correct_links * 1.0 / inferred_links
else:
    ppv = 100.0
print "%.2f" % sens, "%.2f" % ppv
