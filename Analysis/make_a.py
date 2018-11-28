import numpy as np
from scipy import sparse
import sys

ch=sys.argv[1]
res=sys.argv[2]
M=np.loadtxt("M_chr"+ch+"_"+res+"K.txt",dtype='int')
U=np.loadtxt("U_chr"+ch+"_"+res+"K.txt",dtype='int')
Y=np.loadtxt("Y_chr"+ch+"_"+res+"K.txt",dtype='int')

meth=sparse.coo_matrix((M[:,2],(M[:,0],M[:,1])))
unmeth=sparse.coo_matrix((U[:,2],(U[:,0],U[:,1])))
methunmeth=sparse.coo_matrix((Y[:,2],(Y[:,0],Y[:,1])))
meth=meth.tocsr()
unmeth=unmeth.tocsr()
methunmeth=methunmeth.tocsr()

x=np.zeros(meth.shape[1])
for i in range(0, len(x)):
    tot=0
    m=0
    for j in range(0, len(x)):
        tot=tot+unmeth[i,j]+meth[i,j]+methunmeth[i,j]+methunmeth[j,i]
        m=m+meth[i,j]+methunmeth[i,j]

    if tot==0:
        x[i]=0
    else:
        x[i]=float(m)/tot
 
np.savetxt("a_chr"+ch+"_"+res+"K.txt",x,"%0.5f")
