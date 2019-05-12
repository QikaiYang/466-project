import Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import Phylo
import numpy as np

from Bio.Phylo.TreeConstruction import DistanceMatrix

#path = "100M1/100,1000,0.00005,.01,medium_gap_pdf,GTR+second,15,2.0,1/R"


path = "1000M1/1000,1000,.0000082,.005,medium_gap_pdf,GTR+second,35,2.0,1.0/R"

for replicate in range(20):
    print(replicate)
    aln = AlignIO.read(open(path+str(replicate)+"/rose.aln.true.fasta"), "fasta")


    t_mat = []


    base_num = {'A':0, 'T':1, 'C':2, 'G':3, '-':4}

    names = []

    for entry in aln:
        t_counts = np.zeros((5,5))
        t_counts[4][4] = 1
        entry.seq = str(entry.seq).replace('--','')
        seq = entry.seq
        for i in range(1, len(seq)):
            t_counts[base_num[seq[i-1]]][base_num[seq[i]]] += 1

        for i in range(5):
            t_counts[i] /= sum(t_counts[i])

        t_mat.append(np.log(t_counts))
        names.append(entry.id)
    t_mat = np.array(t_mat)
    n = len(aln)
    d_mat1 = []
    d_mat2 = []
    print('calc distance')
    for i in range(n):
        temp1 = []
        temp2 = []
        for j in range(0,i+1):
            '''
            s1 = aln[i].seq
            s2 = aln[j].seq

            s1_prob=0

            for k in range(1,len(s1)):
                s1_prob += t_mat[j][base_num[s1[k-1]]][base_num[s1[k]]]

            s2_prob=0
            for k in range(1, len(s2)):
                s2_prob += t_mat[i][base_num[s2[k-1]]][base_num[s2[k]]]
            distance1 = s1_prob+s2_prob
            '''
            distance2 = np.linalg.norm(t_mat[i]-t_mat[j],ord=2)

            #temp1.append(distance1)
            temp2.append(distance2)
        #d_mat1.append(temp1)
        d_mat2.append(temp2)


    #distance_matrix1 = DistanceMatrix(names, d_mat1)

    #distance_matrix1.format_phylip(open(path+str(replicate)+"/d_matMM.phy", 'wt'))

    distance_matrix2 = DistanceMatrix(names, d_mat2)

    distance_matrix2.format_phylip(open(path+str(replicate)+"/d_matl2.phy", 'wt'))

'''
calculator = DistanceCalculator('ident')
constructor_nj = DistanceTreeConstructor(calculator, 'nj')
tree_nj = constructor_nj.build_tree(aln)

print(tree_nj)
Phylo.write(tree_nj, 'out_tre', 'newick')
'''
