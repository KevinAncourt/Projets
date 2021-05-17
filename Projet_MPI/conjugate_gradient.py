"""
Algorithme du gradient conjugue preconditionne.
"""
import numpy as np
from math import sqrt

def solver_gc( A, b, x0 = None, M = None, tol = 1.E-7, niterMax=100,
               prodMatVect = None, prodPrecVect = None, prodScal = None,
               verbose = False , sommetInt=[], sommetBord={},list_sommet=[]):
    """
    Utilisation :
    importer ce module dans votre ficher python avec :
        from conjugate_gradient import *
    Puis lors de son utilisation :
    sol, info = solver_gc(spMatrix, b, None, tol=1.E-14, M=None, verbose=True)
    print(f"info  : {info}")
    ou sol contient la solution, et info est un tuple contenant le residu final et le nombre d'iteration necessaire pour converger.
    De facon generale, les parametres de la fonction sont :
    A est la matrice du systeme lineaire A.x = b, b le second membre. Le troisieme argument x0 peut etre mis a None si il est au depart nul.
    tol est la convergence desiree, M le preconditionneur, niterMax le nombre max d'iterations tolere, prodMatVect si besoin est une fonction redefinissant le produit
    matrice-vecteur (sinon on prend A.dot comme produit matrice-vecteur), prodPrecVect le produit preconditionneur-vecteur (sinon on utilise M.dot ... ), prodScal si
    on souhaite modifier le produit scalaire, et enfin verbose qui si il est a True, affiche les residus a chaque iteration, et sinon
    n'affiche rien...
    """
    prodA = prodMatVect     
    if prodA is None:
        prodA = A.dot
    prodM = prodPrecVect
    if prodM is None:
        if M is not None:
            prodM = M.dot
    dotS = prodScal
    if dotS is None:
        dotS = np.dot

    if x0 is not None:
        r = b - prodA(A,x0,sommetBord,list_sommet)
    else:
        r = b
    if prodM:
        p = prodM(r)
    else:
        p = r
    z = p

    if x0:
        x = x0.copy()
    else:
        x = np.zeros(b.shape[0], b.dtype)
    err0 = dotS(r,r,sommetInt,sommetBord)
    if verbose :
        print(f"Norme L2 initiale de l'erreur : {sqrt(err0)}")
    err = err0
    nit = 0
    while (err>err0*tol*tol) and (nit < niterMax):
        Apk = prodA(A,p,sommetBord,list_sommet)
        rkzk = dotS(r,z,sommetInt,sommetBord)
        alpha = rkzk/dotS(Apk,p,sommetInt,sommetBord)
        x = x + alpha * p
        r = r - alpha * Apk
        if prodM:
            z = prodM(r)
        else:
            z = r
        beta = dotS(r,z,sommetInt,sommetBord)/rkzk
        p = z + beta * p
        err = dotS(r,r,sommetInt,sommetBord)
        nit += 1
        if verbose:
            print(f"||r_{nit:03}||_L2 : {sqrt(err/err0)}")
    return x, (sqrt(err/err0), nit)