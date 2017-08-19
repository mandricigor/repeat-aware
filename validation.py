


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



#scaffold = [["a", "b"]]
#scaf_orient = [["-", "+"]]

#reference = [["a", "b"]]
#ref_orient = [["+", "+"]]


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
    for scaf in a:
        scaf2 = []
        ori = []
        for contig in scaf:
            contig = contig.split(":::")[:2]
            name = contig[0].split("_")[1]
            scaf2.append(name)
            if contig[1] == "1":
                ori.append("+")
            else:
                ori.append("-")
        answer.append(scaf2)
        orient.append(ori)
    return answer, orient







scaffold, scaf_orient = parse_scaf_file(sys.argv[1])
reference, ref_orient = parse_scaf_file(sys.argv[2])


# paper example
#scaffold = [["f", "a", "b", "c", "d", "e", "d"], ["g"]]
#scaf_orient = [["+", "+", "+", "-", "+", "+", "+"], ["-"]]

#reference = [["a", "b", "c", "b", "d", "e"], ["f", "g"]]
#ref_orient = [["+", "+", "+", "-", "+", "+"], ["+", "+"]]


#print scaffold
#print reference


# dictionary containing relationship: contig -> [reference]
refdict = {}
for count, ref in enumerate(reference):
    for r in ref:
        if r not in refdict:
            refdict[r] = set()
        refdict[r].add(count)
#print refdict

print "OVERALL LINKS:", sum(map(lambda x: len(x) - 1, scaffold))


allscaf = [item for sublist in scaffold for item in sublist]

#print allscaf, "ALLSCAF"

allref = [item for sublist in reference for item in sublist]

#print allref, "ALLREF"

sorient = [item for sublist in scaf_orient for item in sublist]
rorient = [item for sublist in ref_orient for item in sublist]


# scaffmatch

#scaffold = [['39', '41'], ['63', '64', '65', '129', '128', '127', '126', '125', '124', '123', '122', '34', '121', '113', '120', '119', '118', '117', '116', '96', '112', '111', '110', '109', '108', '107', '106', '105', '104', '103', '102', '101', '100', '99', '98', '97', '95', '93', '91', '89', '88', '87', '86', '85', '84', '83', '82', '81', '80', '79', '78', '77', '76', '75', '74', '73', '72', '71', '70', '69', '68', '67', '66', '62', '61', '57', '56', '53', '50', '49', '48', '47', '46', '45', '44', '43', '42', '40', '38', '37', '36', '35', '33', '32', '31', '30', '29', '28', '27', '26', '25', '24', '23', '22', '21', '20', '19', '18', '17', '16', '9', '15', '14', '170', '169', '168', '167', '166', '165', '164', '163', '162', '160', '158', '157', '155', '153', '152', '151', '150', '149', '148', '147', '146', '145', '115', '144'], ['1', '6'], ['156'], ['172'], ['51', '55', '60'], ['114'], ['171', '174', '173'], ['54', '58', '59'], ['4'], ['52'], ['154'], ['2'], ['13', '12', '11', '10', '8', '7', '5'], ['175'], ['161', '159'], ['130', '131', '132', '133', '134', '135', '136', '137', '138', '139'], ['3'], ['92', '90'], ['140', '141', '142', '143']]

#scaf_orient = [['+', '-'], ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '-', '-', '+', '+', '-', '-', '+', '-', '+', '+', '-', '-', '-', '-', '-', '-', '+', '-', '-', '-', '-', '-', '-', '+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-'], ['-', '-'], ['+'], ['+'], ['+', '+', '+'], ['+'], ['+', '-', '+'], ['+', '+', '+'], ['+'], ['+'], ['+'], ['+'], ['+', '-', '-', '-', '+', '+', '+'], ['+'], ['-', '-'], ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-'], ['+'], ['+', '+'], ['-', '-', '-', '-']]


# opera

#scaffold = [['61', '59', '56', '55', '53', '51', '50', '48', '47', '46', '45', '44', '43', '42', '38', '37', '36', '35', '33', '26', '25', '24', '23', '22', '21', '19', '18', '17', '16', '15', '14', '170', '169', '168', '167', '80', '81', '82', '83', '84', '86', '87', '88', '89', '90', '92', '91', '92', '93', '95', '97', '98', '100', '102', '103', '104', '105', '106', '107', '108', '109', '111', '112', '116', '117', '118', '119', '120', '121', '122', '123', '126', '124', '127', '128', '129'], ['40'], ['65'], ['171', '174', '173'], ['49'], ['5'], ['13'], ['1'], ['6'], ['2'], ['58'], ['64', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '166', '165', '164', '163', '162', '160', '157', '155', '153', '152', '151', '150', '149', '148', '147', '146', '145', '144'], ['34'], ['101'], ['31'], ['20'], ['156'], ['30'], ['3'], ['125'], ['4'], ['39'], ['130', '131', '132', '133', '134', '136', '137', '138', '139'], ['113'], ['63'], ['29'], ['115'], ['41'], ['85'], ['52'], ['172'], ['27'], ['54'], ['140', '141', '142', '143'], ['161'], ['60'], ['9'], ['8'], ['32'], ['57'], ['110'], ['135'], ['114'], ['7'], ['159'], ['28'], ['94'], ['175'], ['158'], ['99'], ['96'], ['62'], ['12'], ['11'], ['10']]

#scaf_orient = [['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '-', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '+', '-', '-', '-', '-'], ['+'], ['+'], ['+', '-', '+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['-', '+', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '+', '+', '+', '+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '+', '-'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['-', '-', '-', '-', '-', '-', '-', '-', '-'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['-', '-', '-', '-'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+']]


# opera lg

#scaffold = [['65', '129', '128', '127', '126', '125', '124', '123', '122', '34', '121', '113', '114', '115', '120', '119', '114', '115', '118', '117', '116', '115', '114', '113', '112', '111', '110', '109', '108', '107', '106', '105', '104', '103', '102', '101', '100', '98', '34', '97', '95', '93', '92', '91', '90', '89', '88', '87', '86', '85', '84', '83', '82', '81', '80', '167', '168', '169', '170', '14', '15', '9', '8', '16', '17', '18', '19', '21', '20', '22', '23', '8', '9', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '115', '38', '42', '40', '39', '41', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '58', '55', '56', '59', '61', '62'], ['10', '8', '9', '11'], ['4', '13', '12'], ['62'], ['12'], ['11'], ['65'], ['49'], ['5'], ['13'], ['1'], ['62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '166', '165', '164', '163', '162', '159', '161', '160', '157', '156'], ['6'], ['2'], ['34'], ['31'], ['156'], ['3'], ['125'], ['113'], ['29'], ['115'], ['115', '147', '113', '148', '149', '150', '151', '152', '154', '153', '155', '156'], ['54'], ['60'], ['9'], ['8'], ['57'], ['114'], ['28'], ['94'], ['175'], ['158'], ['62', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '65'], ['99'], ['96'], ['62', '65', '144', '115', '114', '113', '145', '146', '115'], ['62', '140', '141', '142', '143', '65'], ['173', '174', '171', '7'], ['172', '3', '7'], ['7']]
#scaf_orient = [['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '+', '-', '-', '-', '+', '+', '+', '-', '-', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '-', '+', '+', '-', '-', '+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '-', '-', '+', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+'], ['+', '+', '-', '+'], ['+', '+', '-'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+', '+', '+', '+', '+', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '+', '+', '+', '+', '+', '+', '-'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['-', '-', '+', '-', '-', '-', '-', '-', '+', '+', '-', '-'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'], ['+'], ['+'], ['-', '-', '+', '+', '+', '+', '-', '-', '+'], ['-', '-', '-', '-', '-', '-'], ['-', '+', '-', '-'], ['+', '-', '-'], ['+']]

#scaffold = [['96'], ['105'], ['21'], ['126'], ['22'], ['49'], ['65'], ['114'], ['11'], ['19'], ['62'], ['109'], ['16'], ['122'], ['97'], ['8'], ['79'], ['9'], ['64'], ['172'], ['63'], ['113'], ['125'], ['156'], ['20'], ['175'], ['101'], ['173'], ['98', '100'], ['66', '67', '68', '69', '70', '71', '72', '73', '75', '76', '77', '78'], ['171', '174'], ['108', '107', '106', '104', '103', '102'], ['143', '142', '141', '140'], ['17', '18'], ['15', '14', '170', '169', '168', '167', '166', '165', '164', '163', '162', '161', '159', '160', '157', '155', '153', '152', '151', '150', '149', '148', '147', '146', '145', '115', '144'], ['48', '47', '46', '45', '44', '43', '40', '41', '39', '38', '37', '36', '35', '34', '33', '32', '31', '30', '29', '27', '26', '25', '24', '23'], ['139', '138', '137', '136', '135', '134', '133', '132', '131', '130'], ['110', '111', '112', '116', '117', '118', '119', '120', '121'], ['50', '52', '53', '56', '57', '54', '58', '60', '61'], ['95', '93', '92', '90', '91', '90', '89', '88', '87', '86', '85', '84', '83', '82', '81', '80'], ['10', '7', '5', '4', '3', '1', '13', '12'], ['129', '128', '127', '124', '123']]

#scaf_orient = [['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+'], ['+', '+'], ['+', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+'], ['+', '-'], ['+', '+', '+', '+', '-', '-'], ['+', '+', '+', '+'], ['+', '+'], ['-', '-', '-', '-', '-', '-', '+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-'], ['-', '-', '-', '-', '-', '-', '-', '+', '-', '-', '-', '-', '-', '-', '-', '+', '-', '-', '+', '-', '-', '+', '-', '+'], ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'], ['-', '+', '+', '+', '+', '+', '+', '+', '-'], ['+', '+', '+', '+', '+', '+', '+', '+', '+'], ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+'], ['-', '+', '+', '+', '-', '+', '+', '-'], ['+', '+', '+', '+', '+']]



#reference = [['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '9', '8', '11', '9', '8', '12', '13'], ['14', '15', '9', '8', '16', '17', '18', '19', '20', '21', '22', '23', '8', '9', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '62', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '94', '96', '97', '34', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '94', '113', '114', '115', '116', '117', '118', '115', '119', '120', '115', '114', '94', '121', '34', '122', '123', '124', '125', '126', '125', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '62', '140', '141', '142', '143', '62', '144', '115', '114', '113', '94', '145', '146', '96', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '160', '159', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170'], ['171', '172', '173', '174', '175']]


#ref_orient = [['-', '-', '+', '-', '-', '-', '-', '-', '+', '+', '-', '+', '+', '-', '+', '+', '-'], ['+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '+', '+', '-', '+', '+', '-', '-', '+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '+', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '+', '+', '-', '-', '-', '-', '+', '+', '+', '-', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '+', '-', '-', '+', '+', '+', '+'], ['-', '+', '-', '+', '-']]





cons1 = []
cons2 = []
condict = {}

variables = []

#orient_weights = []

i = 0
for ref in reference:
    for r in ref:
        j = 0
        scafn = 0
        c1 = []
        for scaf in scaffold:
            for sc in scaf:
                if r == sc:
                    pos = "%s_%s" % (scafn, j)
                    #print i, pos, r
                    variables.append((i, pos))
                    c1.append((i, pos))
                    p = condict.get(pos, [])
                    p.append(i)
                    condict[pos] = p
                j += 1
            scafn += 1
        #print "-------------"
        i += 1
        cons1.append(c1)
        #orient_weights.extend(o1)

########print orient_weights


#print cons1
#print len(cons1)


#print condict, len(condict)

for x, y in condict.items():
    if len(y) > 1:
        cons = []
        for yy in y:
            cons.append((yy, x))
        cons2.append(cons)

#print cons2

#print len(cons2)





# inter contig constraints
scaf_outer = set()
for scaf in scaffold:
    for a1, a2 in pairwise(scaf):
        scaf_outer.add(tuple(sorted([a1, a2])))


#print scaf_outer



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





#print indices
#print outer_indices, "OUTER INDICES"
#print outer_orient, "OUTER ORIENT" 



cons_pairs_outer = {}
cons_pairs_orient = {}
j = 0
for ref, ori in zip(reference, ref_orient):
    for i in range(len(ref) - 1):
        link = tuple(sorted([ref[i], ref[i + 1]]))
        o1, o2 = ori[i], ori[i + 1]
        if link in scaf_outer:
            #print link, " in scaf outer", j, cons1[j], cons1[j + 1]
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

#print cons_pairs_outer, "CONS PAIRS OUTER"

xz_corr = {} # correspondence between x variables and z
used_x = set()
z = 0
for link in cons_pairs_outer:
    conses = cons_pairs_outer[link]
    conses_orient = cons_pairs_orient[link]
    cind = outer_indices[link]
    cind_orient = outer_orient[link]
    #print link, conses, conses_orient, cind, cind_orient
    zeds = []
    for cons, cori in zip(conses, conses_orient):
        zz = []
        ww = []
        for index, sori in zip(cind, cind_orient):
            if cori == sori:
                w = 1
                ww.append(w)
                cpx.variables.add(lb = [0], ub = [1], types = ["B"], names = ["Z#%s" % z])
                #print cons[0][index[0]], "JJJJJ"
                print cons[0], "CONS0"
                print cons[1], "CONS1"
                print index, "INDEX"
                print link, "LINK"
                inds = ["X#%s#%s" % cons[0][index[0]], "X#%s#%s" % cons[1][index[1]], "Z#%s" % z]
                used_x.add("X#%s#%s" % cons[0][index[0]])
                used_x.add("X#%s#%s" % cons[1][index[1]])
                xz_corr["Z#%s" % z] = ("X#%s#%s" % cons[0][index[0]], "X#%s#%s" % cons[1][index[1]])
                #print inds
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
    #print zeds, "MMMMMMMMMMMMMMMM"
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



#print used_x
#print len(variables), len(used_x)



#print allzeds, zedweights
#print len(allzeds), len(zedweights)


for zed, w in zip(allzeds, zedweights):
    cpx.objective.set_linear(zed, w)

#for var in variables:
    #if var not in used_x:
    #cpx.objective.set_linear("X#%s#%s" % var, 0.0000001)
    #cpx.objective.set_linear("X#%s#%s" % var, 1)


#for var, ori in zip(variables, orient_weights):
#    cpx.objective.set_linear("X#%s#%s" % var, ori)



cpx.objective.set_sense(cpx.objective.sense.maximize)
cpx.set_problem_type(cpx.problem_type.MILP)
#self._cpx.parameters.mip.tolerances.mipgap.set(0.01)
cpx.write("program.txt", filetype="lp")
cpx.solve()
#print cpx.solution.get_values()
#print cpx.variables.get_names()

#for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
#    print x, y

answer_dict = {}
for x in cpx.variables.get_names():
    if "X" in x:
        x = x.split("#")[1]
        answer_dict[x] = []

#for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
#    if y > 0.1 and "X" in x:
#        x = x.split("#")
#        a = x[1]
#        b = x[2].split("_")[0]
#        print "%s %s" % (a, b)



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

print scaf_answer

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

print "AAAAAAAAAAAAAAAAA", scaf_answer


scaf_answer_human = []
asc = iter(allscaf)
for scaf in scaf_answer:
    myscaf = []
    for sc in scaf:
        myscaf.append(asc.next())
    scaf_answer_human.append(myscaf)


print "MMMMMMMMMMMM", scaf_answer_human



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

print corrected_chains

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



for aa, bb, cc, dd in wrong_links:
    is_jumping = False
    is_wrong_ord_ori = False
    is_wrong_copy = False
    is_wrong_ref = False
    print aa, bb, cc, dd
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
        print "AAA - wrong copy"
    elif is_wrong_ref:
        wrong_ref += 1
        print "BBB - wrong ref"
    elif is_jumping and is_wrong_ord_ori:
        jumping_wrong_ord_ori += 1
        print "CCC - jumping + ord_ori"
    elif is_jumping:
        jumping += 1
        print "DDD - jumping"
    elif is_wrong_ord_ori:
        wrong_ord_ori += 1
        print "EEE -wrong ord ori"

print "COPY:" , wrong_copy
print "REF:", wrong_ref
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

#print scaf_answer

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
