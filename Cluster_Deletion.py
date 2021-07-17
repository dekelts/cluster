#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Need to run Cluster_Deletion-prep.py first to generate all cases.

import sys
import scipy.optimize
from copy import deepcopy
from math import ceil

# Compute the branching number for a vector V
def compute(V, s):
	if len(V) <= 1:
		return
	def h(x):
		return sum([x**(-v) for v in V])-1
	X = scipy.optimize.brenth(h,1, 100)
	print s,("%5.3f" % (ceil(X*1000)/1000)),V
	return X

def prod(A, B):
	C = []
	for a in A:
		C += [a+b for b in B]
	return C

# Return a list of adjecent vertices given the adjaceny vector of a vertex
def adj(V):
	return [i for i in range(len(V)) if V[i] == 1]

# Build a graph from a level representation
#
# A level represntation is a tuple (BAmat, BBmat, BB1mat, B1B1mat, B1B2deg)  where
# BAmat[b][a] = 1 if there is an edge between b ∈ B and a ∈ A
# BBmat[b][b2] = 1 if there is an edge between b,b2 ∈ B
# BB1mat[b][c] = 1 if there is an edge between b ∈ B and c ∈ B1
# B1B1mat[c][c2] = 1 if there is an edge between c,c2 ∈ B1
# B1B2deg[c] = |N(c) ∩ B2| for c ∈ B1
#
# The returned representation is a vector G in which G[x] is data of vertex x
# G[x][0] = the level of x (the index i such that x ∈ B_i)
# G[x][1] = the neighbors of x.
#
# The vertices are numbered as follows: First the 3 vertices of A, then the vertices of B,
# the vertices of B1, and the vertices of B2. 
# 
def build_graph(G0):
	BAmat, BBmat, BB1mat, B1B1mat, B1B2deg = G0
	nB = len(BAmat)
	nB1 = len(B1B2deg)
	nB2 = sum(B1B2deg)
	n = 3+nB+nB1+nB2 # |A|=3
	G = [[None,[]] for i in range(n)]

	# Add edges between A and B
	for b in range(nB):
		for a in adj(BAmat[b]):
			bp = 3+b	# The b-th vertex in B is the bp-th vertex in the ordering of all vertices
			G[bp][1].append(a)
			G[a][1].append(bp)

	# Add edges between vertices of B
	for b in range(nB):
		for b2 in adj(BBmat[b]):
			if b < b2:
				bp = 3+b
				b2p = 3+b2
				G[bp][1].append(b2p)
				G[b2p][1].append(bp)

	# Add edges between B and B1
	for b in range(nB):
		for c in adj(BB1mat[b]):
			bp = 3+b
			cp = 3+nB+c
			G[bp][1].append(cp)
			G[cp][1].append(bp)

	# Add edges between vertices of B1
	for c in range(nB1):
		for c2 in adj(B1B1mat[c]):
			if c < c2:
				cp = 3+nB+c
				c2p = 3+nB+c2
				G[cp][1].append(c2p)
				G[c2p][1].append(cp)

	# Add vertices in B1 according to B1B2deg
	dp = 3+nB+nB1	# an index to the current vertex inserted to B2
	for c in range(nB1):
		for k in range(B1B2deg[c]):
			cp = 3+nB+c
			G[cp][1].append(dp)
			G[dp][1].append(cp)
			dp += 1

	# Initialize the level of the vertices of A
	for i in range(3):
		G[i][0] = -1

	# Compute the levels of the vertices of G
	calc_level(G)
	return G

# Compute the levels of the vertices of G (only compute the level if it is <= 1)
def calc_level(G):
	for i in range(3, len(G)):
		G[i][0] = None
	B = []
	for a in range(3):
		for b in G[a][1]:
			if G[b][0] == None:
				G[b][0] = 0
				B.append(b)

	for b in B:
		for c in G[b][1]:
			if G[c][0] == None:
				G[c][0] = 1

# Return the indices of the vertices in level 0 (the set B)
def get_B(G):
	return [i for i in range(len(G)) if G[i][0] == 0]

# Return the indices of the vertices in level 1 (the set B1)
def get_B1(G):
	return [i for i in range(len(G)) if G[i][0] == 1]

# Compute the set F for an edge between a vertex b ∈ B and a vertex c ∈ B1
def calcF(G, b, c):
	if G[b][0] != 0 or G[c][0] != 1:
		print "error!"
	F = []
	# Add edges incident on b
	for x in G[b][1]:
		if G[x][0] == -1:
			F.append((b,x))
		elif x != c and x not in G[c][1]:
			F.append((b,x))
	# Add edges incident on c
	for x in G[c][1]:
		if G[x][0] == None:
			F.append((c,x))
		elif x != b and x not in G[b][1]:
			F.append((c,x))
	return F

# Find an edge between a vertex b ∈ B and a vertex c ∈ B1 with |F_e| >= 3
# If no such edge exists, return None,None
def get_F3_edge(G):
	b2 = c2 = None
	for b in get_B(G):
		for c in G[b][1]:
			if G[c][0] == 1:
				x = len(calcF(G, b, c))
				if x >= 3:
					if x >= 4:
						print "F=4:", x,b,c
					return b,c
#					if b2 == None:
#						b2 = b
#						c2 = c
#
	return b2,c2

def has_B2_neighbors(G, c):
	V = [d for d in G[c][1] if G[d][0] == None]
	return V != []

def find_triangles(G, label = "triangles"):
	for c in get_B1(G):
		if has_B2_neighbors(G,c):
			T = [c2 for c2 in G[c][1] if G[c2][0] == 1 and has_B2_neighbors(G,c2)]
			if T != []:
				print label, len(T), c, ":",
				for c2 in T:
					print c2, len([x for x in G[c2][1] if G[x][0] == 0]),"|",
				print


# def get_F3_edge2(G):
	# for b in get_B(G):
		# for c in G[b][1]:
			# if G[c][0] == 1 and len([x for x in G[c][1] if G[x][0] == None]) >= 2:
				# return b,c
	# return None,None

# delete an edge
def delete_edge(G,b,c):
	del G[b][1][G[b][1].index(c)]
	del G[c][1][G[c][1].index(b)]

# Compute the size of E1
def calc_E1_size(G):
	r = 0
	for c in get_B1(G):
		for x in G[c][1]:
			if G[x][0] == None: # the level of x must be 2 since x is adjacent to c
				r += 1
	return r

# Compute the branching vector for the graph G
# Note: depth is only for debug!
def get_vector(G, depth = 0):
	if verbose or depth == 0:
		for i in range(len(G)):
			print i,G[i]
		print "depth =", depth
		#find_triangles(G, "triangles0")

	b,c = get_F3_edge(G)
	if b == None: # handle the case j > 0
		p = calc_E1_size(G)
		alpha = 1

		u = w = v1 = v2 = False
		for b in get_B(G):
			X = [x for x in G[b][1] if x <= 2]
			if X == [0]: u = True		# u = True iff there is a vertex b∈B s.t. N(b) ∩ A = {u}
			if X == [2]: w = True		# w = True iff there is a vertex b∈B s.t. N(b) ∩ A = {w}
			# if X == [1,2]: v1 = True	# v1 = True iff there is a vertex b∈B s.t. N(b) ∩ A = {v,w}
			# if X == [0,1]: v2 = True	# v2 = True iff there is a vertex b∈B s.t. N(b) ∩ A = {u,v}

		case = 1
		if u and w:
			alpha = 2
			case = 2
		# elif v1 and v2:
			# # There are vertices b,b2 ∈ B s.t. N(b) ∩ A = {v,w} and N(b2) ∩ A = {u,v}
			# # Note that b,b2 can be adjacent or non-adjacent
			# alpha = 1
			# case = 1

		find_triangles(G)
		V = [len(G[c][1]) for c in get_B1(G)]
		#if V != [] and max(V) >= 5:
		print "deg", V

		if verbose: print "p = %d alpha = %d case = %d depth = %d" % (p,alpha,case,depth)
		if verbose: print [alpha+x for x in R[p]]
		return [alpha+x for x in R[p]]

	# V is a vector holding the result
	V = []

	if verbose: print "edge",b,c

	# delete the edge (b,c)
	G2 = deepcopy(G)
	delete_edge(G2,b,c)
	calc_level(G2)
	V += [1+x for x in get_vector(G2, depth+1)]

	# Delete F_{b,c}
	G2 = deepcopy(G)
	F = calcF(G, b, c)
	for x,y in F:
		delete_edge(G2,x,y)
	calc_level(G2)
	if verbose: print "F",b,c,len(F)
	
	V += [len(F)+x for x in get_vector(G2, depth+1)]
	return V

verbose = False

# Compute the vector R
R = [None]*10
R[0] = [0]
R[1] = [1,3]
for i in range(2, 10):
	R[i] = prod(R[i-1],R[1])

G_list = []
if len(sys.argv) > 1:
	verbose = True

if 0:
	BAmat = [[1, 0, 0], [1, 0, 0], [0, 0, 1], [0, 0, 1]]
	BBmat = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
	BB1mat = [[1, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0], [0, 0, 0, 1, 0, 1]]
	B1B1mat = [[0, 1, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1], [0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0]]
	B1B2deg = [0, 0, 1, 0, 0, 1]
	G_list.append(build_graph((BAmat, BBmat, BB1mat, B1B1mat, B1B2deg)))
else:
	lines = file("cd.txt").readlines()
	print len(lines)
	for i in range(0, len(lines), 6):
	#for i in range(0, 30000, 6):
		BAmat = eval(lines[i])
		BBmat = eval(lines[i+1])
		BB1mat = eval(lines[i+2])
		B1B1mat = eval(lines[i+3])
		B1B2deg = eval(lines[i+4])
		G_list.append(build_graph((BAmat, BBmat, BB1mat, B1B1mat, B1B2deg)))

X = [1]*4
for G in G_list:
	print "------------------------------------------------------"
	V = get_vector(G)
	V2 = [None]*4
	V2[0] = V
	V2[1] = [2+x for x in V]+[2+x for x in V]
	V2[2] = [3,6]+[2+x for x in V]
	V2[3] = [5,8,5,8]+[2+x for x in V]
	for i in range(4):
		X[i] = max(X[i],compute(V2[i], "v"+str(i)))

for i in range(4):
	print X[i]
