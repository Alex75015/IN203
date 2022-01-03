# include <chrono>
# include <random>
# include <cstdlib>
# include <sstream>
# include <string>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <mpi.h>
# include <utility>
# include <functional>
#include <typeinfo>
#include <thread>

// Attention , ne marche qu'en C++ 11 ou supérieur :
double approximate_pi( unsigned long nbSamples ) 
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = beginning.time_since_epoch();
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution <double> distribution ( -1.0 ,1.0);
    unsigned long nbDarts = 0;
    // Throw nbSamples darts in the unit square [-1 :1] x [-1 :1]
    for ( unsigned sample = 0 ; sample < nbSamples ; ++ sample ) {
        double x = distribution(generator);
        double y = distribution(generator);
        // Test if the dart is in the unit disk
        if ( x*x+y*y<=1 ) nbDarts ++;
    }
    // Number of nbDarts throwed in the unit disk
    double ratio = double(nbDarts)/double(nbSamples);
    return 4*ratio;
}

template<class Function, class ...Args> // fonction avec pack de paramètres pour fonction à plusieurs var passée en parametre
auto estimate_function_time(Function function, Args... arguments) // ou à la place d'auto : std::pair<std::chrono::high_resolution_clock::duration, long double> mais impose un type long double ce qu'on ne veut pas
{
	auto starting_time = std::chrono::high_resolution_clock::now(); // auto pour type automatique je suppose
	auto result = function(arguments...); // atention au type, ici on a besoin parfois de double ou long double
	auto elasped_time = std::chrono::high_resolution_clock::now() - starting_time;
	return std::make_pair(elasped_time, result);
}

double calcul_pi_non_bloquant_processus_maitre(int tag,int nbp)
{
	double resultat = 0.0;
	double buf=0.0;
	double pi;
	MPI_Request request;
	for(int i =1;i<nbp;i++)
	{
		MPI_Request request;
		// ## attention à bien mettre Irecv et non IRecv !! idem pour Isend et non ISend
		MPI_Irecv(&buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD,&request); // attention le resulat est stocké dans buf
		if (buf == 0.0) // j'étais obligé de faire ça car la première valeur de buf n'étais pas envoyée assez rapidement et donc il ajoutait 0 à résultat au lieu de 3.1xx
		{
			MPI_Wait(&request,MPI_STATUS_IGNORE);
			resultat = resultat + buf;
		}
		else 
			resultat = resultat + buf;
		
		//~ std::cout << "resultat " << resultat << "\n";
		MPI_Wait(&request,MPI_STATUS_IGNORE); //besoin de & ou non?
	}
	
	//~ std::cout << "resultat " << resultat << "et 13xpi = " << 3.14*13 << " \n";
	
	pi = resultat/(nbp-1); // attention le processus 0 ne calcule pas
	
	std::cout << "L'approximation du nombre pi avec réception non bloquante vaut : " << pi << ".\n";
	return pi;
}

int calcul_pi_non_bloquant(int nargs, char* argv[], unsigned long nbSamples)
{

	//~ MPI_Init( &nargs, &argv );

	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);

	int nbp;
	MPI_Comm_size(globComm, &nbp);
	
	//~ std::cout << "nbre " << nbp << "\n";

	int rank;
	MPI_Comm_rank(globComm, &rank);

	int tag = 78;
	double pi;
	
	 if (rank == 0)
    {
		auto pi = estimate_function_time(calcul_pi_non_bloquant_processus_maitre,tag, nbp);
		std::cout << "L'approximation de pi avec réception non bloquante a été fait en " << pi.first.count() << " ticks.\n";

	}
	
	else
	{
		double buf = approximate_pi(nbSamples);
		//~ std::cout << "L'approximation local du nombre pi vaut : " << buf << " pour le processus " << rank << ".\n";
		MPI_Status status;
		MPI_Send(&buf, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	
	
	//~ MPI_Finalize();
	return EXIT_SUCCESS;
}

double calcul_pi_bloquant_processus_maitre(int tag,int nbp)
{
	double resultat = 0.0;
	double buf=0.0;
	double pi;
	MPI_Status status;
	for(int i =1;i<nbp;i++)
	{
		double buf;
		MPI_Recv(&buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status); // attention le resulat est stocké dans buf
		resultat = resultat + buf;
	}
	
	//~ std::cout << "resultat " << resultat << "et 13xpi = " << 3.14*13 << " \n";
	
	pi = resultat/(nbp-1); // attention le processus 0 ne calcule pas
	
	std::cout << "L'approximation du nombre pi avec récéption bloquante vaut : " << pi << ".\n";
	return pi;
}

int calcul_pi_bloquant(int nargs, char* argv[], unsigned long nbSamples)
{

	//~ MPI_Init( &nargs, &argv );

	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);

	int nbp;
	MPI_Comm_size(globComm, &nbp);
	
	//~ std::cout << "nbre " << nbp << "\n";

	int rank;
	MPI_Comm_rank(globComm, &rank);

	int tag = 78;
	//~ double pi;
	
	 if (rank == 0)
    {
		
		//~ pi = calcul_pi_non_bloquant_processus_maitre(tag, nbp);
		auto pi = estimate_function_time(calcul_pi_bloquant_processus_maitre,tag, nbp);
		std::cout << "L'approximation de pi avec réception bloquante a été fait en " << pi.first.count() << " ticks.\n";

	}
	
	else
	{
		double buf = approximate_pi(nbSamples);
		//~ std::cout << "L'approximation local du nombre pi vaut : " << buf << " pour le processus " << rank << ".\n";
		MPI_Status status;
		MPI_Send(&buf, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	
	
	//~ MPI_Finalize();
	
	return EXIT_SUCCESS;
}

int main( int nargs, char* argv[] )
{
	MPI_Init( &nargs, &argv );

    
    calcul_pi_non_bloquant(nargs, argv, 10000);
    
	calcul_pi_bloquant(nargs, argv, 10000);
	
	MPI_Finalize();
	
	return 0;
}


//## A faire :
//## utiliser le code du TP9 de IN204 pour comparer le temps de calcul_pi avec et sans envoie bloquant
