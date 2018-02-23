from simLib import simulateTree, simulateSNPs
from random import shuffle

def main():
    n = 40
    N = 10**6
    alpha = 0 #float(7)/4/N
    mu = float(40)/4/N
    myTree = simulateTree(n,N,alpha)
    M = simulateSNPs(myTree,mu)
    myTree.write_to_path("Simulated_Tree.nwk","newick")

    # randomly swap the order of the mutations
    # i.e. the rows of M and output 
    #the transposed version of M (which is the proper form of a SNP matrix)  

    Si = [0]*n
    for c in M:
        Si[sum(c)-1] += 1

    with open("Xi.txt",'w') as f:
        for i in range(n):
            f.write(str(i+1) + " " + str(Si[i])+ "\n")

    columns = list(range(len(M)))
    shuffle(columns)

    with open("SNP_Matrix.txt","w") as f:
        for i in range(n):
            for j in columns:
                f.write(str(M[j][i]))
            f.write("\n")    

if __name__=='__main__':
    main()        
