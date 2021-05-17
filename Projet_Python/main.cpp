#include "maillage.hpp"
#include <time.h>

// Question 1:

void test_affichage_maillage(const maillage & Th){
    cout << Th << endl; // affiche l'ensemble du maillage
}

// Question 2:

void test_triangle_adjacent(maillage & Th){
  // Boucle affichant le tableau contenant l'ensemble des triangles adjacents pour chaque triangle du maillage Th
  
    int * tadj = Th.triangles_adj(Th);
    for (int i = 0; i < 3*Th.nt; i += 3) {
        cout << " Le triangle " <<i/3 << " a pour triangles adjacents : ";
        for (int j = 0 ; j < 3 ; j++ ) {
            cout << tadj[i + j] << ' ';
         }
         cout << endl;
      }
    delete tadj;
}

void test_affichage_triangle_adjacent(maillage & Th,int triangle_depart,int sommet){
    // Renvoie le triangle adjacent du triangle "triangle_depart" opposé au sommet "sommet"
  
    int * tadj = Th.triangles_adj(Th);
    cout << "Le triangle adjacent au triangle "<< triangle_depart << " opposé au sommet "<< sommet << " est : " <<Th.return_triangle_adj(tadj, triangle_depart, sommet) << endl;
    delete tadj;
}

//Question 3:

void test_algorithme1(maillage & Th,R2 p){
    // Algorithme recherchant le triangle du maillage Th qui contient le point p
  
    int k = Th.Algorithme_1(0,p);
    cout << " Le point (" << p << ") appartient au triangle "<< k << endl;
}

//Question 4 :

void test_algorithme1_generaliser(maillage & Th,maillage & Th2){
    // Affiche l'ensemble des sommets du maillage Th2 ainsi que le triangle du maillage Th contenant chaque sommet
    // Affiche un message si le sommet n'appartient à aucun triangle
  
    int *  tab = Th.Q4(Th2);
    for (int i=0;i<Th2.nv; i++){
        if (tab[i]== -1){
            cout << "Le sommet " << i << " n'est pas inclus dans le maillage"<<endl;
        }
        else{
            cout << "le sommet "<< i << " est inclus dans le triangle "<<tab[i] <<endl;
        }
    }
    delete tab;
}

int main(int argc,const char ** argv)
{
    assert( argc>2); // Vérifie qu'on envoie deux arguments à l'exécution
    maillage Th(argv[1]);
    maillage Th2(argv[2]);
    
    R2 p(-0.6,-0.3); // définit un point de R2
    
    // Mesure du temps
    clock_t start, end;
    double t;
    srand(time(NULL));
    start = clock();

    //Ensemble des tests:
    
    //test_affichage_maillage(Th);
    //test_triangle_adjacent(Th);
    //test_affichage_triangle_adjacent(Th,1,2);
    //test_algorithme1(Th,p);
    test_algorithme1_generaliser(Th,Th2);
    
    end = clock();
    t = ((double)end - start) / CLOCKS_PER_SEC;
    cout << " temps en secondes : " << t << endl;

    
    
    return 0;
}
