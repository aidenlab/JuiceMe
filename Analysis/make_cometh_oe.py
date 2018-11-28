import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

import sys

# chromosome, resolution are passed in
ch=sys.argv[1]
res=sys.argv[2]

# load meth-meth, unmeth-unmeth, meth-unmeth matrices
M=np.loadtxt("M_chr"+ch+"_"+res+"K.txt",dtype='int')
U=np.loadtxt("U_chr"+ch+"_"+res+"K.txt",dtype='int')
Y=np.loadtxt("Y_chr"+ch+"_"+res+"K.txt",dtype='int')
# load WGBS probability vector
a=np.loadtxt("a_chr"+ch+"_"+res+"K.txt")

E=np.zeros(shape=(len(a),len(a)))
O=np.zeros(shape=(len(a),len(a)))
Q=np.zeros(shape=(len(a),len(a)))

# create sparse matrices
meth=sparse.coo_matrix((M[:,2],(M[:,0],M[:,1])))
unmeth=sparse.coo_matrix((U[:,2],(U[:,0],U[:,1])))
unmethmeth=sparse.coo_matrix((Y[:,2],(Y[:,0],Y[:,1])))
meth=meth.tocsr()
unmeth=unmeth.tocsr()
unmethmeth=unmethmeth.tocsr()
for i in range(0, len(a)):
    for j in range(0, len(a)):
        total=unmeth[i,j]+meth[i,j]+unmethmeth[i,j]+unmethmeth[j,i]
        E[i,j]=((a[i]*a[j])+((1-a[i])*(1-a[j])))*total
        O[i,j]=meth[i,j]+unmeth[i,j]
        Q[i,j]=total

np.savetxt("observed_chr"+ch+"_"+res+"K.txt", O, "%0.5f")
np.savetxt("expected_chr"+ch+"_"+res+"K.txt", E, "%0.5f")
np.savetxt("total_chr"+ch+"_"+res+"K.txt", Q, "%0.5f")
