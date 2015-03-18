#!/usr/bin/python

# Sept 21, 2012
# Michael Soltys
# Converting graph G into shuffle string w_G

import sys
import math
import itertools
import re

# When matrix is provided in a file matrix.txt,
# read in the matrix from the file & print to standard output

def ReadMatrixFromInput(matrixname):
	mymatrix = open(matrixname,'r')
	Q=[0]
	R=[]
	for line in mymatrix.readlines():
		print line[:-1]
		Q = Q+[int(line[i]) for i in range(len(line)-1)]
		R = R+[Q]
		Q = [0]
	R = [[0 for i in range(len(line))]]+R
	mymatrix.close
	return R

# If e=(v_i,v_j) is an edge, check that 
# v_i^1,v_i^2,v_j^1,v_j^2 have the proper values

def CheckConnected(V,i,j):
	if (separate == 0):
		if (V[i][1]<V[j][1] and V[j][1]<V[i][2] and V[i][2]<V[j][2]):
			return 1
		elif (V[j][1]<V[i][1] and V[i][1]<V[j][2] and V[j][2]<V[i][2]):
			return 1
		else:
			return 0

	#
	# when separate == 1 we allow this situation; i.e., we do not insist
	# on overlap
	#      __          __
	#     /  \        /  \
	#    /    \      /    \
	# --+------+----+------+--
	# v_i^1 v_i^2 v+j^1  v_j^2
	#
	elif (separate == 1):
		if (V[i][1]<V[j][1] and V[i][2]<V[j][2]):
			return 1
		elif (V[j][1]<V[i][1] and V[j][2]<V[i][2]):
			return 1
		else:
			return 0

# If e=(v_i,v_j) is NOT an edge, check that 
# v_i^1,v_i^2,v_j^1,v_j^2 have the proper values

def CheckNotConnected(V,i,j):

	# first the strict mode, laxness == 0, where the arches are not
	# allowed to overlap at the end points; that is, given two points
	# that are not connected in the graph, their corresponding arches
	# are such that one of them is properly contained in the other:
	#        _____
	#       / ___ \
	#      / /   \ \
	#     / /     \ \
	#    + +       + +
	#
	# note that if -per is used, then the search is strict a priori; thus,
	# it is a waste to do "-per -lax" together, as the 'laxness' is not
	# used in a "-per" search.
	if (laxness == 0):
		if (V[i][1]<V[j][1] and V[j][1]<V[j][2] and V[j][2]<V[i][2]):
			return 1
		elif (V[j][1]<V[i][1] and V[i][1]<V[i][2] and V[i][2]<V[j][2]):
			return 1
		else:
			return 0

	# then we have the lax mode, -lax, laxness == 1, where the arches
	# are allowed to meet at the end points, so for example we can have:
	#        _____
	#       /_____\              __     __
	#      /       \            /  \   /  \
	#     /         \          /    \ /    \
	#    +           +   or   +      +      +  etc...
	#
	elif (laxness == 1):
		if (V[i][1]<=V[j][1] and V[j][1]<=V[j][2] and V[j][2]<=V[i][2]):
			return 1
		elif (V[j][1]<=V[i][1] and V[i][1]<=V[i][2] and V[i][2]<=V[j][2]):
			return 1
		elif (V[i][1]<V[i][2] and V[i][2]==V[j][1] and V[j][1]<V[j][2]):
			return 1
		elif (V[j][1]<V[j][2] and V[j][2]==V[i][1] and V[i][1]<V[i][2]):
			return 1
		else:
			return 0

# Given a list of integers e, construct from them a candiate
# assignment V

def ConstructCandidate(e):
	
	# This is a nice example of list merge (first operator '+')
	# and list comprehension - where we construct the [0,e[i],e[i+1]]
	# it is here that one can see the inventor of Python to be
	# a mathematician - comprehension is a nice concept from logic

	V = [[0,0,0]]+[[0,e[i],e[i+1]] for i in range(0,P-2,2)]
	return V

# Verify that the given candidate, V, works for the matrix M, where N
# is the number of vertices

def VerifyCandidate(M,N,V):
	
	for i in range(1,N):
		for j in range(i+1,N+1):
			if (M[i][j]==0 and CheckNotConnected(V,i,j)==0):
				return 0
			elif (M[i][j]==1 and CheckConnected(V,i,j)==0):
				return 0
	return 1

# While searching for a good V, this gauge gives the progress
# in percentage of iterations/permutations checked.

def CheckGauge(p,total,search_type):

	accurate = float(p)/float(total)

	# if we are examining iterations
	if (search_type == 0):
		sys.stdout.write("Computing %d iterations: %d%c \r" % 
		(total,int(round(accurate,2)*100),"%"))
		
	# if we are examining permutations
	if (search_type == 1):
		sys.stdout.write("Computing %d permutations: %d%c \r" % 
		(total,int(round(accurate,2)*100),"%"))

# Build a matrix of zeros, of size (N+1)x(N+1),
# Our matrices will then be embedded in the NxN 
# principal submatrix

def ConstructZeroMatrix(N):

			M = []
			for i in range(N+1):
				row = [0 for j in range(N+1)]
				M = M+[row]

			return M

# After the search for V is completed, construct corresponding
# string, or inform that it cannot be done.

def DisplayOutput(found,V,N):

	if (found == 1):
		print
		print "Found good string:"
		for i in range(1,N+1):
			w[V[i][1]]=w[V[i][2]]=S[i]
			print "%d-%d" % (V[i][1],V[i][2])
		print ''.join(w)
	else:
		print
		print "Could not find good string."

# The main search subroutine

def Search(M,N,P,gauge,search_type):

	# now we define the search space

	if (search_type == 0): # plain iterations
		
		# Cartesian product: {1,...,P-1}^(P-1)
		allvalues = [range(1,P) for i in range(1,P)]

		# counter of the number of iterations examined
		iters_total = math.pow(P-1,P-1)

		SearchSpace = itertools.product(*allvalues)

	elif (search_type == 1): # permutations

		# counter of the number of permutations examined
		iters_total = math.factorial(P-1)

		SearchSpace = itertools.permutations(range(1,P))
	
	# 'found' is flag for having found the right V; initially, we don't
	# have it, so we haven't found it. 'iters_counter' counts the numbe
	# of iterations of the loop - used with CheckGauge
	found = 0
	iters_counter = 0

	for e in SearchSpace:
		
		# check progress if "-g" switch was set
		iters_counter += 1
		if (gauge == 1):
			CheckGauge(iters_counter,iters_total,search_type)

		# get the candidate corresponding to permutation e
		V=ConstructCandidate(e)

		# verify the candidate
		if (VerifyCandidate(M,N,V) == 1):
			found = 1
			break

	DisplayOutput(found,V,N)

def ConstructMatrix(N,b):
	k = 0
	for i in range(1,N):
		for j in range(i+1,N+1):
			M[i][j] = b[k]
			k += 1

	return M

def PrintMatrix(M,N):

	for i in range(1,N+1):
		print M[i][1:]

def PrintMatrixWithDegree(M,N,deg):

	for i in range(1,N+1):
		print M[i][1:],deg[i]

def ComputeDegrees(M,N):
	deg = [0 for j in range(0,N+1)]

	# In order to compute the degree of the i-th vertex -- keeping in mind
	# that (our undirected) graphs are represented as upper-triangular
	# 0-1 matrices with 0s on the main diagonal -- we add up the
	# following entries:
	#
	#      i-th column
	# +-----------------------+
	# |    x                  |
	# |    x                  |
	# |    x                  |
	# |     xxxxxxxxxxxxxxxxxx| i-th row
	# |                       |
	# .                       .

	for i in range(1,N+1):
		deg_value = 0
		
		# sum up elements in i-th column, start at row 1, end at row i-1
		# the (i,i) entry is always zero so ignore it
		for j in range(1,i):
			deg_value += M[j][i]

		# sum up elements in i-th row, start at column i+1, end at column
		# N; again, the (i,i) entry is always zero so ignore it
		for k in range(i+1,N+1):
			deg_value += M[i][k]

		# enter value of degree of i-th vertex, and move on to i+1 vertex,
		# by first zeroing deg_valu out
		deg[i] = deg_value

	return deg

def CheckDegreeCond(deg,N):
	for i in range(1,N):
		if (deg[i]<deg[i+1]):
			return 0

	return 1
	

## FIRST WE DEFINE ALL THE GLOBAL PARAMETERS
## these are usually flags for different search options

# this string variable contains the alphabet
# first entry is empty to start counting at 1
# this implicitly limits size of matrix to 26x26 -
# this is fine as program slows down a lot at 6x6 already
S = " abcdefghijklmnopqrstuvwxyz"

## We accept the following switches:
# 	-g 		gauge 
#		-per		permutations only search type
#		-lax	laxness
#		-sep	allow separate arches for connected vertices

## We give default values to the global variables associated with the
## switches

gauge = 0	# '-g' switch turns 'gauge = 1'

# for search type:
# nothing - means plain iterations:	search_type = 0
# -per means permutations							search_type = 1
search_type = 0

# laxness refers to types of V we consider
# nothing - means strict mode:				laxness = 0
# -lax means V can meet at end points	laxness = 1
laxness = 0	# '-lax' switches to 'laxness = 1'

# if separate == 1 it means that we do not insist that v_j^1<v_i^2 or
# that v_i^1<v_j^2; that is, we allow "separate arches" in the case
# that (v_i,v_j) are connected. That is, we allow the following
# situation:
#      __          __
#     /  \        /  \
#    /    \      /    \
# --+------+----+------+--
# v_i^1 v_i^2 v+j^1  v_j^2
#
# default is separate == 0, i.e., we insist on overlapping arches for
# connected vertices.
separate = 0	# '-sep' switch changes to 'separate == 1'

# this is only meaningful in the batch mode, i.e., when we are
# visiting -r i-j. By default all graphs ixi...jxj are examined, but
# the switch '-deg' turns 'degree == 1' and only some of the graphs in
# those ranges are examined, namely all those graphs with the "degree
# property": a graph G has this property if for all the nodes in
# V={1,2,3,...,k} (where we examine the Gs such that i<=k<=j) we have
# that deg(1)>=deg(2)>=deg(3)>=...>=deg(k). The idea is that the
# remaining graphs, those that do not have the degree property, are
# isomorphic to some graph of same size V which does have the degree
# property.
degree = 0 		# '-deg' switch changes to 'degree == 1'

## First, we examine the command line, in order to parse
## the input and classify the switches
## we also start the program - all auxiliary function definitions have
## been given above

# number of arguments from the command line; remember:
# sys.argv[0] = name of the program ('shuffle.py' in our case)
# sys.argv[1] = main switch: -m input matrix; -r for a range 
# sys.argv[2] = name of the file containing the matrix or range
# sys.argc[3] = 1st switch
# sys.argc[4] = 2nd switch ... etc

number_arguments = len(sys.argv)

if (number_arguments >= 3): 

	# Check if there are more than 2 arguments, i.e., check
	# if there are auxiliary switches
	if (number_arguments >= 4):
		for i in range(3,number_arguments):
			if (sys.argv[i] == "-g"):
				gauge = 1
			if (sys.argv[i] == "-per"):
				search_type = 1
			if (sys.argv[i] == "-lax"):
				laxness = 1
			if (sys.argv[i] == "-sep"):
				separate = 1
			if (sys.argv[i] == "-deg"):
				degree = 1

	if (sys.argv[1] == "-m"):
		# The input is a particular matrix given in 3rd argument
    # M is the matrix that describes the graph
		M = ReadMatrixFromInput(sys.argv[2])

		# Once we have M, we can initialize some global variables
	
		# remember M is padded with a row & column of zeros
		# this padding makes vertex 1 to be given by row 1
		# (otherwise, without padding, vertex 1 would be given
		# by row 0)
		# thus, N = number of vertices, is given as:
		N = len(M)-1 
	
		# P gives the length of the resulting string
		# we increase by 1 so that for convenience we number
		# the positions of the string starting at 1
		P = 2*N+1

		# this list variable holds P spots that will be populated
		# with the alphabet symbols on the way to create final string
		w = ['' for i in range(P)]

		# start search for string corresponding to M
		Search(M,N,P,gauge,search_type)

	if (sys.argv[1] == "-r"):
		# Search over all matrices of size nxn to mxm, inclusive, where
		# the 3rd argument is n-m
		match = re.match('([0-9]*)-([0-9]*)',sys.argv[2])
		low_range = int(match.group(1))
		high_range = int(match.group(2))

		# number of matrices examined:
		matrix_counter = 1

		for N in range(low_range,high_range+1):
			P = 2*N+1

			M = ConstructZeroMatrix(N)

			allvalues = [range(2) for i in range((N*N-N)/2)]

			for b in itertools.product(*allvalues):
				if (degree == 1):

					M = ConstructMatrix(N,b)
					
					deg = ComputeDegrees(M,N)

					if (CheckDegreeCond(deg,N) == 1):
						print "Matrix #: %d" % (matrix_counter)
						
						PrintMatrixWithDegree(M,N,deg)
						
						matrix_counter += 1

						w = ['' for i in range(P)]
						Search(M,N,P,gauge,search_type)
						print "------------------------------------"	
					
				else:
					print "Matrix #: %d" % (matrix_counter)
					matrix_counter += 1
					
					M = ConstructMatrix(N,b)
				
					PrintMatrix(M,N)

					w = ['' for i in range(P)]
					Search(M,N,P,gauge,search_type)
					print "------------------------------------"	

else:
	print "shuffle.py is a Python program"
	print "Michael Soltys - September 21, 2012"
	print "RCS revision 1.7"
	print "------------------------------------------------------------"
	print "Description:"
	print "------------------------------------------------------------"
	print "shuffle.py implements a function that attempts to compute"
	print "a reduction from clique to shuffle. It does so as follows:"
	print "on input an adjacency matrix of an undirected graph G, i.e.,"
	print "a matrix over Z_2^{NxN}, it computes a string w_G, with the"
	print "property that any vertex v is mapped to a pair of indices"
	print "(v1,v2) corresponding to the same symbol in w_G. Note that"
	print "given G=(V,E), if |V|=N, then |w_G|=2N, and w_G is a string"
	print "over an alphabet of N distinct symbols, each occurring"
	print "exactly twice. Further, (u,v) in E iff the arches that"
	print "correspond to (u1,u2) and (v1,v2) overlap in a legal way."
	print "The result is that the string w_G represents the connections"
	print "of G, and it is such that it contains a shuffle of size k"
	print "iff the graph has a clique of size k. Note however, that"
	print "the reduction does not always work."
	print "------------------------------------------------------------"
	print "Usage:"
	print "------------------------------------------------------------"
	print "There are two main ways of using this program:"
	print "\t shuffle.py -m matrix.txt"
	print "\t shuffle.py -r 2-5"
	print "The first example extracts a Z_2^{NxN} matrix from file"
	print "matrix.txt and searches for the corresponding string"
	print "The second example computes the strings of all matrices"
	print "over Z_2 of sizes 2x2,3x3,4x4,5x5. Of course, we can use"
	print "-r n-m for any n <= m, but keep in mind that m=6 and"
	print "bigger slows down noticeably"
	print "------------------------------------------------------------"
	print "Switches:"
	print "------------------------------------------------------------"
	print "The main operation has to go first; then there are a number"
	print "of switches that can be employed:"
	print "\t -g \t shows a gauge of progress"
	print "\t -per \t searches for a solution over permutations only"
	print "\t -lax \t allows for lax mode; without it is strict mode"
	print "\t -sep \t allows for separate arches for connected points"
	print "\t -deg \t only search graphs s.t., d(1)>=d(2)>=...>=d(k)"
	print "------------------------------------------------------------"
	print "More detailed explanation of the switches:"
	print "------------------------------------------------------------"
	print "-lax:"
	print "We have a strict (by default, without -lax) and lax"
	print "(with -lax) mode when checking (u,v) not in E."
	print "In the strict mode the arches are not allowed"
	print "to overlap at the end points; that is, given 2 points that are"
	print "not connected in the graph, their corresponding arches are such"
	print "that one of them is properly contained in the other:"
	print "       _____     "
	print "      / ___ \    "
	print "     / /   \ \   "
	print "    / /     \ \  "
	print "   + +       + + "
	print 
	print "With the -lax switch the arches are allowed to meet at the end" 
	print "points, so for example we can have:"
	print "       _____                                    "
	print "      /_____\              __     __            "
	print "     /       \            /  \   /  \           "
	print "    /         \          /    \ /    \          "
	print "   +           +   or   +      +      +  etc... "
	print
	print "Note that with -per, the search is strict a priori; thus, it is"
	print "it is a waste to do -per -lax together, as the 'laxness' is not"
	print "used in a -per search."
	print "------------------------------------------------------------"
	print "-sep:"
	print "Switch -sep means that we do not insist that v_j^1<v_i^2 or"
	print "that v_i^1<v_j^2; that is, we allow separate arches in the case"
	print "that (v_i,v_j) are connected. That is, we allow the following"
	print "situation:"
	print "     __          __      "
	print "    /  \        /  \     "
	print "   /    \      /    \    "
	print "--+------+----+------+-- "
	print "v_i^1 v_i^2 v+j^1  v_j^2 "
	print 
	print "the default is to insist on overlapping arches for"
	print "connected vertices. The difference is that without -sep,"
	print "the solution to a 5x5 clique is:"
	print "1-6"
	print "2-7"
	print "3-8"
	print "4-9"
	print "5-10"
	print "abcdeabcde"
	print "while with -sep the solution to a 5x5 clique is:"
	print "1-2"
	print "3-4"
	print "5-6"
	print "7-8"
	print "9-10"
	print "aabbccddee"
	
	sys.exit()
