# !/usr/bin/env

from dendropy import Tree,Node,Taxon
from random import random, randint,shuffle, expovariate
from math import exp
from numpy.random import poisson, exponential
from sys import stdout

def simulateTreeTopology(n):
# simulate a binary tree of n leaves
    leaves = [Node()]
    nodeOrder = []
    myTree = Tree(seed_node = leaves[0])

    for i in range(n-1):
        r = randint(0,i)
        a = Node()
        b = Node()
        p = leaves[r]
        p.add_child(a)
        p.add_child(b)
        leaves[r] = a
        leaves.append(b)
        nodeOrder.append(p)

    IDs = list(range(1,n+1))
    i = 0
    shuffle(IDs)   
    for leaf in leaves:
        leaf.taxon = Taxon(label=str(IDs[i]))
        leaf.time = 0
        i += 1        

    return myTree,nodeOrder    

def coalescentRate(N,n):
    return float(n*(n-1)/2)/N

def shrinkPopulation(N,alpha):
    return N*exp(-alpha)


def simulateTreeBranches(tree,nodeOrder,N,alpha):
# simulate branches following coalescent model
# with current population size N and exponential growth parameterized by alpha
    currTime = 0
    currN = 2*N

    for k in range(len(nodeOrder)+1,1,-1):    
        currNode = nodeOrder[k-2]

        if alpha == 0:
            currTime += exponential(1/coalescentRate(currN,k))
        else:    
            while 1:
                r = random()
                p = coalescentRate(currN,k)
                if r < p: # successed
                    break
                currTime += 1
                currN = shrinkPopulation(currN,alpha)    

        currNode.time = currTime
        for child in currNode.child_node_iter():
            child.edge_length = currTime - child.time     

def simulateTree(n,N,alpha):
    tree,nodeOrder = simulateTreeTopology(n)
    simulateTreeBranches(tree,nodeOrder,N,alpha)
    return tree

def orderTreeNodes(myTree):
# randomly choosing a valid ordering of the nodes in myTree
# this is helpful when we want to fix the tree topology but 
# shuffle the coalescent/speciation events 

    activeNodes = [myTree.seed_node]
    nodeOrder = [myTree.seed_node]

    while activeNodes:
        r = randint(0,len(activeNodes)-1)
        p = activeNodes[r]
        nodeOrder.append(p)

        # remove p from activeNodes by swapping it to the end and pop out
        q = activeNodes[-1]
        activeNodes[r] = q
        activeNodes.pop()

        for c in p.child_node_iter():
            if not c.is_leaf():
                activeNodes.append(c)
            else:
                c.time = 0 
                       
    return nodeOrder    

def simulateTreeFromTopology(myTree,N,alpha):
    nodeOrder = orderTreeNodes(myTree)
    simulateTreeBranches(myTree,nodeOrder,N,alpha)

def simulateSNPs(tree,mu):
# "drop" mutations onto the simulated tree
# and output the SNP matrix
# note: the output matrix M here is a transposed version
# of a standard SNP matrix (for convenient implementation)
    M = []
    n = len(list(tree.leaf_node_iter()))
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.list = [node.taxon.label]
        else:
            node.list = []
            for c in node.child_node_iter():
                node.list += c.list
                c.list = None
        if node.edge_length is None:
            break    
        k = poisson(mu*node.edge_length)
        if k <= 0:
            continue
        column = [0]*n
        for x in node.list:
            column[int(x)-1] = 1
              
        for i in range(k):
            M.append(column)

    return M

