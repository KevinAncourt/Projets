#!/bin/env python
# -*- coding: utf-8 -*-
import mesh
import fem
import laplacian
import splitter
from math import cos,sin,pi,sqrt
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import VisuSplitMesh as VSM
import VisuSolution as VS
from scipy.sparse import linalg as sp_linalg
from mpi4py import MPI
from conjugate_gradient import *



comm= MPI.COMM_WORLD
rank= comm.rank
size=comm.size

def prodMatVect(mat,vect,bord,list_sommet):
    y=np.zeros(vect.shape[0])
    
    for i in range(len(list_sommet)):
        for j in range(len(list_sommet)):
            y[i] +=mat[i,j] * vect[j]
       
    for i in bord:
        voisin=[]
        for j in bord[i]:
            voisin.append(y[j])
        comm.send(voisin,dest=i)
        
    for i in bord:
        value_voisin=comm.recv(source =i)
        for elt_recv in bord[i]:
            y[elt_recv] += value_voisin[bord[i].index(elt_recv)] 
    return y
    



def prodScal(x,y,sommet_interne,bord):
    scal=0.
    list_bord=[]
    for i in bord:
        for j in bord[i]:
            if j not in list_bord:
                list_bord.append(j)
    
    for i in sommet_interne:
        scal += x[i] * y[i]
    
    for j in list_bord:

        cout=1
        for k in bord:
            for l in bord[k]:
                if l ==j:
                    cout +=1

        scal += (1./cout)*(x[j] * y[j])

    scalglobal = comm.allreduce(scal, MPI.SUM)
    return scalglobal
 
    

def g(x,y) :
    return cos(2*pi*x)+sin(2*pi*y)


m = mesh.read("CarrePetit.msh")
coords    = m[0]
elt2verts = m[1]
nbVerts = coords.shape[0]
nbElts  = elt2verts.shape[0]
begVert2Elts, vert2elts = mesh.compvert2elts(elt2verts)


begRows, indCols = fem.comp_skel_csr_mat(elt2verts, (begVert2Elts, vert2elts) )
nz = begRows[-1]

spCoefs = np.zeros( (nz,), np.double)

for iElt in range(nbElts):
    iVertices = elt2verts[iElt,:]
    crd1 = coords[iVertices[0],:]
    crd2 = coords[iVertices[1],:]
    crd3 = coords[iVertices[2],:]
    matElem = laplacian.comp_eltmat((crd1[0],crd1[1]), (crd2[0],crd2[1]),
                                    (crd3[0],crd3[1]))
    fem.add_elt_mat_to_csr_mat((begRows,indCols,spCoefs),
                               (iVertices, iVertices, matElem))

# Assemblage second membre :
f = np.zeros(nbVerts, np.double)
for iVert in range(nbVerts):
    if ( coords[iVert,3] > 0 ) :
        f[iVert] += g(coords[iVert,0],coords[iVert,1])
b = np.zeros(nbVerts, np.double)
for i in range(nbVerts) :
    for ptR in range(begRows[i],begRows[i+1]):
        b[i] -= spCoefs[ptR]*f[indCols[ptR]]
        
# Il faut maintenant tenir compte des conditions limites :
for iVert in range(nbVerts):
    if coords[iVert,3] > 0: # C'est une condition limite !
        # Suppression de la ligne avec 1 sur la diagonale :
        for i in range(begRows[iVert],begRows[iVert+1]):
            if indCols[i] != iVert :
                spCoefs[i] = 0.
            else :
                spCoefs[i] = 1.
        # Suppression des coefficients se trouvant sur la colonne iVert :
        for iRow in range(nbVerts):
            if iRow != iVert :
                for ptCol in range(begRows[iRow],begRows[iRow+1]):
                    if indCols[ptCol] == iVert :
                        spCoefs[ptCol] = 0.
                        
        b[iVert] = f[iVert]
# On definit ensuite la matrice :
spMatrix = sparse.csc_matrix((spCoefs, indCols, begRows),
                             shape=(nbVerts,nbVerts))


# Résolution du problème séquentiel avec gradient conjugue :
iteration = 1
def display_iter(x_k):
    global iteration
    iteration += 1

# Résoud à l'aide du gradient conjugué.
sol,info = sp_linalg.cg(spMatrix, b, x0=None, tol=1.E-14, M=None,callback=display_iter)
d = spMatrix.dot(sol) - b




nbDoms = size
ndsDomains = splitter.node( nbDoms, coords )

i = 0
for a in ndsDomains :
    i += 1

ndsDomains = splitter.node( nbDoms, coords )
ndsDomains_Proc=ndsDomains[rank]
nb_sommet=len(ndsDomains_Proc)

etsDomains = splitter.element( nbDoms, (elt2verts, coords) )
etsDomains_Proc=etsDomains[rank]
nb_element=len(etsDomains_Proc)


########################################################
################# Début de la parallèlisation #################

tri_global=[]
list_sommet=[]
Coords=[]


for e in etsDomains_Proc:
	tri_global.append(list(elt2verts[e])) #elt2vert version local
	list_sommet.append(elt2verts[e][0])
	list_sommet.append(elt2verts[e][1])
	list_sommet.append(elt2verts[e][2])

	
list_sommet=list(set(list_sommet)) #liste des sommets sans doublons et dans l'ordre croissant
list_sommet.sort()
#print(f" Pour le proc  {rank} list sommet est  {list_sommet}")

	
for i in range (len(list_sommet)):
	Coords.append((list(coords[list_sommet[i]]))) # Liste des coordonnées locales 


for i in range(len(tri_global)):
	for j in range(3):
		tri_global[i][j]= list_sommet.index(tri_global[i][j])  # Numérotation locale


##############################################################
######################  Création interface ########################

sommets=comm.allgather(sendobj=list_sommet) # on reuni la liste des sommets 
bord={i:[] for i in range(size) if i!= rank} 
sommet_interne = [list_sommet.index(list_sommet[i]) for i in range(len(list_sommet)) ]

for i in range(len(list_sommet)):
    for j in range(size):
        if j != rank:
            if list_sommet[i] in sommets[j]:
                bord[j].append(list_sommet.index(list_sommet[i]))
                if list_sommet.index(list_sommet[i]) in sommet_interne:
                    sommet_interne.remove(list_sommet.index(list_sommet[i]))
                

elt2doms = np.zeros((nb_element,), np.double)
ia = 0.	
for a in range(len(etsDomains_Proc)):
	elt2doms[a] = rank
    

# Calcul l'interface :
ie = 0
mask = np.array([-1,]*nbVerts, np.short)
for e in tri_global :
    d = elt2doms[ie]
    if mask[e[0]] == -1 :
        mask[e[0]] = d
    elif mask[e[0]] != d :
        mask[e[0]] = -2
    if mask[e[1]] == -1 :
        mask[e[1]] = d
    elif mask[e[1]] != d :
        mask[e[1]] = -2
    if mask[e[2]] == -1 :
        mask[e[2]] = d
    elif mask[e[2]] != d :
        mask[e[2]] = -2
    ie += 1

nbInterf = 0
for m in mask :
    if m == -2 :
        nbInterf += 1

interfNodes = np.empty(nbInterf, np.long)
nbInterf = 0
for im in range(mask.shape[0]):
    if mask[im] == -2 :
        interfNodes[nbInterf] = im
        nbInterf += 1
        
  
#Conversion liste en array        
Coords =np.array(Coords)
tri_global= np.array(tri_global)
        

##############################################################
######################  Création matrice et vecteur locales ########################
begVert2Elts_new, vert2elts_new = mesh.compvert2elts(tri_global)
begRows_new, indCols_new = fem.comp_skel_csr_mat(tri_global, (begVert2Elts_new, vert2elts_new) )
nz_new = begRows_new[-1]

spCoefs_new = np.zeros( (nz_new,), np.double)

for iElt in range(nb_element):
    iVertices = tri_global[iElt,:]
    crd1 = Coords[iVertices[0],:]
    crd2 = Coords[iVertices[1],:]
    crd3 = Coords[iVertices[2],:]
    matElem_new = laplacian.comp_eltmat((crd1[0],crd1[1]), (crd2[0],crd2[1]),
                                    (crd3[0],crd3[1]))
    fem.add_elt_mat_to_csr_mat((begRows_new,indCols_new,spCoefs_new),
                               (iVertices, iVertices, matElem_new))
                               
nbVerts_new = len(list_sommet)
                               
# Assemblage second membre :
f = np.zeros(nbVerts_new, np.double)
for iVert in range(nbVerts_new):
    if ( Coords[iVert,3] > 0 ) :
        f[iVert] += g(Coords[iVert,0],Coords[iVert,1])
b = np.zeros(nbVerts_new, np.double)

for i in range(nbVerts_new) :
    for ptR in range(begRows_new[i],begRows_new[i+1]):
        b[i] -= spCoefs_new[ptR]*f[indCols_new[ptR]]

# Il faut maintenant tenir compte des conditions limites :
for iVert in range(nbVerts_new):
    if Coords[iVert,3] > 0: # C'est une condition limite !
        # Suppression de la ligne avec 1 sur la diagonale :
        for i in range(begRows_new[iVert],begRows_new[iVert+1]):
            if indCols_new[i] != iVert :
                spCoefs_new[i] = 0.
            else :
                spCoefs_new[i] = 1.
        # Suppression des coefficients se trouvant sur la colonne iVert :
        for iRow in range(nbVerts_new):
            if iRow != iVert :
                for ptCol in range(begRows_new[iRow],begRows_new[iRow+1]):
                    if indCols_new[ptCol] == iVert :
                        spCoefs_new[ptCol] = 0.
                   
        b[iVert] = f[iVert]

##############################################################
#  Echange à l'interface pour que le second membre soit un restriction du vecteur globale #
for i in bord:
    list_b=[]
    for j in bord[i]:
        list_b.append(b[j])
    comm.send(list_b,dest=i)
        
for i in bord:
    value_b=comm.recv(source =i)
    for elt_recv in bord[i]:
        b[elt_recv] += value_b[bord[i].index(elt_recv)] 

# On definit ensuite la matrice :
spMatrix_new = sparse.csc_matrix((spCoefs_new, indCols_new, begRows_new),
                             shape=(nbVerts_new,nbVerts_new))


sol, info = solver_gc(spMatrix_new, b, None, tol=1.E-14, M=None, verbose=True,prodScal = prodScal, prodMatVect=prodMatVect ,sommetInt=sommet_interne, sommetBord=bord,list_sommet=list_sommet)

print(f"||r_final|| = {info[0]}, nombre itérations : {info[1]}")
#print(sol)

d = spMatrix_new.dot(sol) - b
print("rank : ", rank ," ||A.x-b||/||b|| = {}".format(sqrt(d.dot(d)/b.dot(b))))


#Visualisation de la solution :
if rank==3:
    VS.view( Coords, tri_global, sol, title = "Solution" )
    
solution=[  0.9999999999999996 , 0.9999999999999996 , 0.9999999999999997 , 0.9999999999999997 , 0.3090169943775894 , -0.8090169943711615 , -0.8090169943788371 , 0.30901699437179997 , 0.04894348370484613 , 0.4122147477004772 , 1.5877852522840126 , 1.9510565162967695 , 0.30901699437494623 , -0.8090169943698257 , -0.8090169943810945 , 0.3090169943699743 , 0.04894348370484613 , 0.4122147477004772 , 1.5877852522840126 , 1.9510565162967695 , 0.24115933881193144 , 0.38914643451072145 , 0.11973709884040552 , 0.18378802620488616 , 0.37282426165065036 , -0.2141514788039678 , 0.5177491177198342 , 0.5830980129442604 , -0.06979682388052677 , 0.12734241721085837 , 0.8665053834998058 , 0.15582737534357413 , 0.7716129591792209 , 0.3329558390258152 , 0.2396742972816799 , 0.16305159139043118 , 0.21762580141883375 , 1.03513538200728 , -0.1581474615917546 , 0.27719080656383194 , -0.1727431909200857 , 1.0096692825554254 , 0.31514909686663856 , 0.028091388931591035 , 0.08100627512540258]
  
sol_tab = comm.gather(sol, 0) 

if rank == 0 :
    sol_refaite = np.zeros(nbVerts)
    for i in range(size) :
        for j in range(len(sol_tab[i])) :
            if sol_refaite[sommets[i][j]] == 0 :
                sol_refaite[sommets[i][j]] = sol_tab[i][j]

    err = np.linalg.norm(sol_refaite-solution, 2)/np.linalg.norm(solution, 2)
    print("erreur en norme l2 ", err)
    
    
