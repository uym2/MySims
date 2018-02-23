#! /usr/bin/env python

from simLib import simulateTreeFromTopology
from sys import argv
from dendropy import Tree

infile = argv[1]
outfile = argv[2]

N = 10**6
alpha = 0

myTree = Tree.get_from_path(infile,'newick')

simulateTreeFromTopology(myTree,N,alpha)

myTree.write_to_path(outfile,'newick')
