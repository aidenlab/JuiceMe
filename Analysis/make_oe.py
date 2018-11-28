import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import sys

# chromosome, resolution are passed in
ch=sys.argv[1]
res=sys.argv[2]

# load meth-meth, unmeth-unmeth, unmeth-meth matrices
M=np.loadtxt("M_chr"+ch+"_"+res+"K.txt",dtype='int')
U=np.loadtxt("U_chr"+ch+"_"+res+"K.txt",dtype='int')
Y=np.loadtxt("Y_chr"+ch+"_"+res+"K.txt",dtype='int')
# load WGBS probability vector
a=np.loadtxt("a_chr"+ch+"_"+res+"K.txt")

# create sparse matrices
meth=sparse.coo_matrix((M[:,2],(M[:,0],M[:,1])))
unmeth=sparse.coo_matrix((U[:,2],(U[:,0],U[:,1])))
unmethmeth=sparse.coo_matrix((Y[:,2],(Y[:,0],Y[:,1])))
meth=meth.tocsr()
unmeth=unmeth.tocsr()
unmethmeth=unmethmeth.tocsr()

# this is unnecessary, should be able to do as sparse matrices
# but this makes it explicit
mm=np.zeros(shape=(len(a),len(a)))
uu=np.zeros(shape=(len(a),len(a)))
mu=np.zeros(shape=(len(a),len(a)))
um=np.zeros(shape=(len(a),len(a)))
for i in range(0, len(a)):
      for j in range(0,len(a)):
              mm[i,j]=meth[i,j]
              uu[i,j]=unmeth[i,j]
              um[i,j]=unmethmeth[i,j]
              mu[i,j]=unmethmeth[j,i]

# correlation based on methylation status
Mcorr=np.zeros(shape=(len(a),len(a)))
Ucorr=np.zeros(shape=(len(a),len(a)))
for i in range(0, len(a)):
    for j in range(0, len(a)):
          if (mm[i,j]+um[i,j])==0:
                Mcorr[i,j]=float('NaN')
          else:
                Mcorr[i,j]=mm[i,j]/(mm[i,j]+um[i,j])
          if (mu[i,j]+uu[i,j])==0:
                Ucorr[i,j]=float('NaN')
          else:
                Ucorr[i,j]=mu[i,j]/(mu[i,j]+uu[i,j])

# average correlation status given locus i
avg_Mcorr=np.zeros(len(a))
avg_Ucorr=np.zeros(len(a))
count_col1=0
count_col2=0
for i in range(0, len(a)):
    sum_col1=0
    sum_col2=0
    for j in range(0,len(a)):
          if not np.isnan(Mcorr[j,i]):
                sum_col1=sum_col1+Mcorr[j,i]
          if not np.isnan(Ucorr[j,i]):
                sum_col2=sum_col2+Ucorr[j,i]
    avg_Mcorr[i]=sum_col1
    avg_Ucorr[i]=sum_col2
    if sum_col1>0:
        count_col1=count_col1+1
    if sum_col2>0:
        count_col2=count_col2+1

# this code removes centromeres/telomeres
tmp=np.delete(Mcorr,np.where(avg_Mcorr==0),0)
tmp2=np.delete(tmp,np.where(avg_Mcorr==0),1)
Mcorr=tmp2
tmp=np.delete(Ucorr,np.where(avg_Ucorr==0),0)
tmp2=np.delete(tmp,np.where(avg_Ucorr==0),1)
Ucorr=tmp2

# plot correlation methylated-unmethylated
f=plt.figure()
fig=plt.imshow(Mcorr-Ucorr, cmap='seismic',interpolation='nearest',vmin=-1,vmax=1)
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.colorbar()
#plt.show()
f.savefig("Mcorr_minus_Ucorr_chr"+ch+"_"+res+"K.pdf")

# positional observed and expected
O=np.zeros(shape=(len(a),len(a)))
E=np.zeros(shape=(len(a),len(a)))
Q=np.zeros(shape=(len(a),len(a)))
for i in range(0, len(a)):
    for j in range(0, len(a)):
        total=uu[i,j]+mm[i,j]+mu[i,j]+um[i,j]
        E[i,j]=a[i] # row vector of WGBS methylation
        # probability that we are methylated given in contact with j
        O[i,j]=mm[i,j]+mu[i,j]
        Q[i,j]=total
        if total==0:
              O[i,j]=float('NaN')
        else:
              O[i,j]=O[i,j]/total

# this code removes centromeres/telomeres
tmp=np.delete(O,np.where(avg_Mcorr==0),0)
tmp2=np.delete(tmp,np.where(avg_Mcorr==0),1)
O=tmp2
tmp=np.delete(E,np.where(avg_Ucorr==0),0)
tmp2=np.delete(tmp,np.where(avg_Ucorr==0),1)
E=tmp2

# plot positional O - E
f=plt.figure(dpi=1000)
fig=plt.imshow(O-E, cmap='seismic',vmin=-0.4,vmax=0.4, interpolation='nearest')
plt.colorbar()
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
#plt.show()
f.savefig("Positional_O_minus_E_chr"+ch+"_"+res+"K.pdf")

# plot positional O / E
f=plt.figure()
fig=plt.imshow(O/E, cmap='seismic', vmax=2)
plt.colorbar()
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
#plt.show()
f.savefig("Positional_O_over_E_chr"+ch+"_"+res+"K.pdf")

# plot average correlation
avg_Mcorr=avg_Mcorr[avg_Mcorr>0]
avg_Ucorr=avg_Ucorr[avg_Ucorr>0]
avg_Mcorr=avg_Mcorr/count_col1
avg_Ucorr=avg_Ucorr/count_col2

f=plt.figure()
p1 = plt.bar(range(count_col2), avg_Ucorr, 0.1, color='r', edgecolor='r')
p2 = plt.bar(range(count_col1), avg_Mcorr-avg_Ucorr, 0.1, color='b', edgecolor='b', bottom=avg_Ucorr)
y=np.zeros(count_col1)
y.fill(np.mean(a))
plt.plot(range(count_col1),y, '-g')

#plt.show()
f.savefig("avg_Mcorr_avg_Ucorr_chr"+ch+"_"+res+"K.pdf")
