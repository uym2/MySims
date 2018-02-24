#! /usr/bin/env python

from simLib import simulateTreeFromTopology
from dendropy import Tree
import argparse
from sys import stdout

parser = argparse.ArgumentParser()

parser.add_argument("-i","--inputFile",required=True,help="The input tree in Newick format")
parser.add_argument("-N","--population",required=True,help="Population size")
parser.add_argument("-g","--growth",required=False,help="Growing rate (exponential) of the population. Default: 0")
parser.add_argument("-o","--outputFile",required=False,help="The name of the output tree. Default: stdout")

args = vars(parser.parse_args())

infile = args["inputFile"]
outfile = args["outputFile"] if args["outputFile"] else None

N = float(args["population"])
alpha = int(args["growth"]) if args["growth"] else 0

myTree = Tree.get_from_path(infile,'newick')

simulateTreeFromTopology(myTree,N,alpha)

if outfile is not None:
    myTree.write_to_path(outfile,'newick')
else:
    stdout.write(myTree.as_string('newick'))    
