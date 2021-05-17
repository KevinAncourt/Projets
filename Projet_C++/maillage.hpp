#include "R2.hpp"
#include <cassert>
#include <vector>
#include <fstream>
#include <random>

using namespace std;

mt19937 gen(109);

class maillage {

public:

  int nv; // Nombre de sommets
  int nt; // Nombre de triangles
    
  vector< vector<double> > v; // Coordonnées des sommets
  vector< vector<int> > t; // Les sommets par triangle
    
  maillage(const char *filename) {
    ifstream f(filename);
    assert(f); // Teste que le fichier est bien ouvert
        
    int unused,label;
    double x,y;
        
    f >> nv >> nt >> unused;
    v.resize(nv);
    t.resize(nt);
       
    assert(f.good());
        
    for(int i=0;i<nv;i++){
      f >> x >> y >> label;
      v[i].push_back(x);
      v[i].push_back(y);
      v[i].push_back(label);
      assert(f.good());
    }
        
    for(int k=0;k<nt;++k){
      int v1,v2,v3,useless;
      f >> v1 >> v2 >> v3 >> useless;
      t[k].push_back(v1-1);
      t[k].push_back(v2-1);
      t[k].push_back(v3-1);
      t[k].push_back(useless);
            
    }
    cout<< " End Read " << nv << " " << nt << endl;
        
  }
    
  int operator()(int k,int i) const {
    // Surcharge de l'opérateur ()
    // Vérifie qu'on ne demande pas un sommet absurde ou un triangle absurde
    assert(k<nt);
    assert(k>=0);
    assert(i>=0);
    assert(i<3);
    return t[k][i];}

  friend ostream& operator << (ostream& o, const maillage &Th);

      
  int * triangles_adj(const maillage & Th) {
    int * t_adj = new int [Th.nt*3];
    int m = Th.nv;
    int nn = Th.nt*3; // majorant du nombre d'arêtes
    int ne = 0;
    std::vector<int> head(m,-1),next(nn),v(nn),ta(nn);
    for(int t=0 ;t<Th.nt;t++){
      for(int e=0 ;e<3 ;e++)
	{
	  t_adj[3*t+e]=-1; // Initialise à -1 pour traiter le cas où on ne trouve pas de triangle adjacent opposé au sommet

	  int i = Th(t, (e+2)%3); // On récupère les sommets opposés au sommet considéré
	  int j = Th(t, (e+1)%3);
        
	  if( j<i) std::swap(i,j);
        
	  int ke=ne;
	  
	  bool exist = false; // permet de ne pas considérer une arête déjà traitée
	  
	  for (int p= head[i]; p>=0; p= next[p]){ // boucle sur les éléments d'une classe d'équivalence
	    if(v[p] == j ) {
	      ke = p;
	      exist=true; // Si l'arête a déjà été traitée on ne la traite pas à nouveau
	      break;}
	  }
	  if( !exist) // Si l'arête n'a pas encore été traitée
	    {
	      next[ne]= head[i];
	      head[i]= ne;
	      v[ne]=j;
	      ta[ne] = 3*t+e; // tableau de tous les triangles
	      ne++;
	    }
	  if( exist )
	    {
	      int tt=ta[ke], ee= tt%3 ;
	      tt /= 3;
	      t_adj[3*t+e]=tt;
	      t_adj[3*tt+ee]=t;
	    }
	}
    }

    return t_adj;
  }

  int return_triangle_adj( int * t_adj, int num_t, int i ) {
    // Renvoie le triangle adjacent de num_t opposé au sommet i
    assert(num_t<nt);
    assert(num_t>=0);
    assert(i>=0);
    assert(i<3);
    return t_adj[ 3 * num_t + i];
  }
    
  int Algorithme_1 (int k, const R2 & p){
      
    int * tadj = (*this).triangles_adj((*this));
    float mes;
    vector<int> triangle_aire_neg(1) ;
    int g = k;

    while (triangle_aire_neg.size() != 0) {
      triangle_aire_neg.resize(0);
      
      for (int i=0; i<3; ++i){
    
          R2 A(v[t[g][i]][0],v[t[g][i]][1]), B(v[t[g][(i+1)%3]][0],v[t[g][(i+1)%3]][1]);
    
          mes = det(A,B,p); // Calcule l'aire du triangle ABp
    
          if ( mes < 0 ){ // On stocke quand les aires sont négatives
              triangle_aire_neg.push_back(return_triangle_adj(tadj,g,(i+2)%3)); }
    }
      
    if(triangle_aire_neg.size()==1){
        g=triangle_aire_neg[0];
        if (g==-1) { // Le cas où le sommet n'appartient à aucun triangle
            break; }
    }
    else if(triangle_aire_neg.size()==2){
        std::bernoulli_distribution R(0.5);
        g=triangle_aire_neg[R(gen)];
        if (g==-1) { // Le cas où le sommet n'appartient à aucun triangle
            break; }
          }
      }
    
    delete tadj;
    return g; }
    
    
  int * Q4 ( const maillage & Th2) {
    // Pour chaque point du maillage Th2, l'algorithme trouve le triangle du maillage Th qui le contient
    
    int *tab = new int [Th2.nv];
    
    for (int i=0; i<Th2.nv; i++){
      
      R2 p(Th2.v[i][0],Th2.v[i][1]);
      tab[i]=(*this).Algorithme_1(1,p);
    }

    return tab;
  }
  
        
};
    

inline ostream& operator << (ostream& o, const maillage &Th){
  
  // Affiche les coordonnées de chaque sommet
  for(int i=0; i<Th.nv; i++){
      o << "sommet numéro " << i << " de coordonnées : ";
    for(int j=0; j<2; j++){
      o << Th.v[i][j] << " ";
    }
    o << endl;

  }
  o << endl;
  
  // Affiche les triangles
  for(int i=0; i<Th.nt; i++){
    o << "triangle numéro " << i << " de sommets : " ;
    for(int j=0; j<3; j++){
      o <<  Th.t[i][j] << " ";
    }
    o << endl;

  }
  o << endl;
  return o;
}







    
