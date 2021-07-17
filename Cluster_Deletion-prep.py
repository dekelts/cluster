#!/usr/bin/pypy
# -*- coding: utf-8 -*-

# Return a list of adjecent vertices given the adjaceny vector of a vertex
def adj(V):
	return [i for i in range(len(V)) if V[i] == 1]

# Return the set {0,...,n-1}\{ex}
def range_exclude(n, ex):
	return [i for i in range(n) if i != ex]

# Return the set {0,...,n-1}\{ex,ex2}
def range_exclude2(n, ex1, ex2):
	return [i for i in range(n) if i != ex1 and i != ex2]

# Return a list of all monotone sequences of length at most n over the alphabet {0,...,sigma-1}
def all_monotone_sequences(n, sigma):
	if n == 1:
		return [[c] for c in range(sigma)]
	A = all_monotone_sequences(n-1, sigma)
	B = []
	for x in A:
		if len(x) == n-1:
			B += [[c]+x for c in range(x[0]+1)]
	return A+B

def all_sequences(n, sigma):
	if n == 0:
		return [[]]
	B = []
	for x in all_sequences(n-1, sigma):
		B += [x+[c] for c in range(sigma)]
	return B

# Generate a list of BAmat vectors, where # BAmat[b][a] = 1 if there is an edge between b ∈ B and a ∈ A
def generate_BAmat_list():
	R = []
	for X in all_monotone_sequences(4, 6):
		# The sequence X is over alphabet of size 6 since a vertex in B has 1 or 2 neighbors in A, so there are 6 cases.
		# We go over monotone sequences to eliminate symetric cases
		if len([c for c in X if c <= 3]) > 2:	# if |Bu| > 2
				continue
		if len([c for c in X if c >= 2]) > 2:	# if |Bw| > 2
				continue
		if len(X) < 0:
			continue

		BAmat = [[0]*3 for x in X]
		for i,c in enumerate(X):
			if c == 0:
				BAmat[i][0] = 1
			elif c == 1:
				BAmat[i][1] = 1
				BAmat[i][2] = 1
			elif c == 2:
				BAmat[i][1] = 1
			elif c == 3:
				BAmat[i][0] = 1
				BAmat[i][2] = 1
			elif c == 4:
				BAmat[i][2] = 1
			elif c == 5:
				BAmat[i][0] = 1
				BAmat[i][1] = 1
		R.append(BAmat)
	return R

# Generate a list of XXmat vectors, where # XXmat[b][b2] = 1 if there is an edge between x,x2 ∈ X
# n is the size of X
def generate_XXmat_list(n):
	R = []
	for X in all_sequences(n*(n-1)/2, 2):
		XXmat = [[0]*n for x in range(n)]
		k = 0
		for i in range(1,n):
			for j in range(i):
				if X[k] == 1:
					XXmat[i][j] = 1
					XXmat[j][i] = 1
				k += 1
		R.append(XXmat)
	return R

# Generate a list of BB1deg vectors, where # BB1deg[b] = |N(b) ∩ B1| for b ∈ B
def generate_BB1deg_list(k):
#	return [[x+1 for x in X] for X in all_sequences(k,2)]
#	The line above assume we delete vertices with empty Nnext
	return [X for X in all_sequences(k, 3) if sum(X) > 0]
#	# |N(b) ∩ B1| ∈ {0,1,2} for every b ∈ B so the alphabet size is 3.
#	# We need at least one index in which X[b] > 0

# An item in the BB1a list is a list of the B1-endpoints of the edges between B and B1 (k = num of edges)
# The edges are ordered according to their B-endpoints.
# For each edge, its B1-endpoint is either the endpoint of a previous edge, or a new vertex.
def generate_BB1a_list(k):
	if k == 0:
		return [[]]
	if k == 1:
		return [[1]]
	r = []
	for X in generate_BB1a_list(k-1):
		m = max(X)
		for y in range(1,(m+1)+1):
			r.append(X+[y])
	return r

# Generate a dictionary: BB1mat_dict[BB1deg] is a list of BB1mat vectors corresponding to BB1deg,
# where BB1mat[b][c] = 1 if there is an edge between b ∈ B and c ∈ B1
def generate_BB1mat_dict():
	BB1a_list = [generate_BB1a_list(k) for k in range(8+1)]
	BB1mat_dict = {}
	for Bsize in range(1,4+1):
		for BB1deg in BB1deg_list[Bsize]:
			key = tuple(BB1deg)
			BB1mat_dict[key] = []
			n = sum(BB1deg) # number of B-B1 edges
			for X in BB1a_list[n]:
				m = max(X) # size of B1
				BB1mat = [[0]*m for i in range(4)]
				k = 0
				for b in range(Bsize):
					for i in range(BB1deg[b]):
						BB1mat[b][X[k]-1] = 1
						k += 1
				if sum([sum(x) for x in BB1mat]) == n:
					BB1mat_dict[key].append(BB1mat)
	return BB1mat_dict

def generate_B1B2deg_list(k):
#	return [[x+1 for x in X] for X in all_sequences(k,2)]
#	The line above assume we delete vertices with empty Nnext
	return all_sequences(k, 3)

#	R = []
#	for X in list(all_strings(k, 3)):
#		if 2 in X:
#			R.append(X)
#	return R

def test1(BAmat, BBmat, BB1deg):
	nB = len(BAmat)
	# check A-B edges
	for b in range(nB):
		for x in adj(BAmat[b]):
			f = 0

			# P3s for the form A-x-b 
			for x2 in range_exclude(3, x):
				xx2edge = 1
				if (x,x2) in [(0,2),(2,0)]:
					xx2edge = 0
				if xx2edge != BAmat[b][x2]:
					f += 1
#			f += 1

			# P3s of the form x-b-B 
			for b2 in range_exclude(nB, b):
				if BBmat[b][b2] != BAmat[b2][x]: f += 1
			# P3s of the form x-b-B1
			f += BB1deg[b]
			if f > 3:
				return False

	# check B-B edges
	for b in range(nB):
		for b2 in adj(BBmat[b]):
			f = 0
			# P3s of the form b-b2-A
			for i in range(3):
				if BAmat[b][i] != BAmat[b2][i]: f += 1
			# P3s of the form b-b2-B
			for b3 in range_exclude2(nB, b, b2):
				if BBmat[b][b3] != BBmat[b2][b3]: f += 1
			# P3s of the form b-b2-B1 (lower bound)
			f += abs(BB1deg[b]-BB1deg[b2])
			if f > 3:
				return False
	return True

def generate_D1_list():
	D1_list = []
	for BAmat in BAmat_list:
		for BBmat in XXmat_list[len(BAmat)]:
			for BB1deg in BB1deg_list[len(BAmat)]:
				if test1(BAmat, BBmat, BB1deg):
					D1_list.append([BAmat, BBmat, BB1deg])
	return D1_list

def test2(BAmat, BBmat, BB1mat):
	nB = len(BAmat)
	nB1 = len(BB1mat[0])

	# check B-B edges
	for b in range(nB):
		for b2 in adj(BBmat[b]):
			f = 0
			# P3s of the form b-b2-A
			for i in range(3):
				if BAmat[b][i] != BAmat[b2][i]:	f += 1
			# P3s of the form b-b2-B
			for b3 in range_exclude2(nB, b, b2):
				if BBmat[b][b3] != BBmat[b2][b3]: f += 1
			# P3s of the form b-b2-B1
			for c in range(nB1):
				if BB1mat[b][c] != BB1mat[b2][c]: f += 1
			if f > 3:
				return False
	return True

def generate_D2_list(D1_list):
	D2_list = []
	for BAmat,BBmat,BB1deg in D1_list:
		for BB1mat in BB1mat_dict[tuple(BB1deg)]:
			if test2(BAmat, BBmat, BB1mat):
				D2_list.append([BAmat, BBmat, BB1mat]) #, B1Bmat])
	return D2_list

def test3(BAmat, BBmat, BB1mat, B1B2deg):
	nB = len(BAmat)

	# check B-B1 edges
	for b in range(nB):
		for c in adj(BB1mat[b]):
			f = 0
			# P3s of the form A-b-c
			f += sum(BAmat[b])
			# P3s of the form b-c-B
			for b2 in range_exclude(nB, b):
				if BBmat[b][b2] != BB1mat[b2][c]: f += 1
			# Note that we don't know the number of P3s of the form b-c-B1
			#
			# P3s of the form A-b-c-B2
			f += B1B2deg[c]
			if f > 3:
				return False
	return True

def generate_D3_list(D2_list):
	D3_list = []
	for BAmat,BBmat,BB1mat in D2_list:
		n = len(BB1mat[0])
		for B1B2deg in B1B2deg_list[n]:
			if test3(BAmat, BBmat, BB1mat, B1B2deg):
				D3_list.append([BAmat, BBmat, BB1mat, B1B2deg])
	return D3_list

def test4(BAmat, BBmat, BB1mat, B1B1mat, B1B2deg):
	nB = len(BAmat)
	nB1 = len(BB1mat[0])

	# check B-B1 edges
	for b in range(nB):
		for c in adj(BB1mat[b]):
			f = 0
			# P3s of the form A-b-c
			f += sum(BAmat[b])
			# P3s of the form b-c-B
			for b2 in range_exclude(nB, b):
				if BBmat[b][b2] != BB1mat[b2][c]: f += 1
			# P3s of the form b-c-B1
			for c2 in range_exclude(nB1, c):
				if BB1mat[b][c2] != B1B1mat[c][c2]: f+= 1
			# P3s of the form A-b-c-B2
			f += B1B2deg[c]
			if f > 3:
				return False

#	return True
	# check B1-B1 edges
	for c in range(nB1):
		for c2 in adj(B1B1mat[c]):
			f = 0
			# P3s of the form c-c2-B
			for b in range(nB):
				if BB1mat[b][c] != BB1mat[b][c2]: f += 1
			# P3s of the form c-c2-B1
			for c3 in range_exclude2(nB1, c, c2):
				if B1B1mat[c][c3] != B1B1mat[c2][c3]: f += 1
			# P3s of the form c-c2-B2 (lower bound)
			f += abs(B1B2deg[c]-B1B2deg[c2])
			if f > 3:
				return False

	return True

def generate_D4_list(D3_list):
	D4_list = []
	for BAmat,BBmat,BB1mat,B1B2deg in D3_list:
		n = len(BB1mat[0]) # size of B1
		for B1B1mat in XXmat_list[n]:
			if test4(BAmat, BBmat, BB1mat, B1B1mat, B1B2deg):
				D4_list.append([BAmat, BBmat, BB1mat, B1B1mat, B1B2deg])
	return D4_list

BAmat_list = generate_BAmat_list()
XXmat_list = [generate_XXmat_list(k) for k in range(6+1)] # |B1| <= 6
BB1deg_list = [generate_BB1deg_list(k) for k in range(4+1)]
BB1mat_dict = generate_BB1mat_dict()
B1B2deg_list = [generate_B1B2deg_list(k) for k in range(8+1)]

D1_list = generate_D1_list()
print "D1:", len(D1_list)
D2_list = generate_D2_list(D1_list)
print "D2:", len(D2_list)
D3_list = generate_D3_list(D2_list)
print "D3:", len(D3_list)
D4_list = generate_D4_list(D3_list)
print "D4:", len(D4_list)

fh  = open("cd.txt", "w")
for BAmat, BBmat, BB1mat, B1B1mat, B1B2deg in D4_list:
	fh.write(str(BAmat)+"\n")
	fh.write(str(BBmat)+"\n")
	fh.write(str(BB1mat)+"\n")
	fh.write(str(B1B1mat)+"\n")
	fh.write(str(B1B2deg)+"\n\n")
fh.close()
