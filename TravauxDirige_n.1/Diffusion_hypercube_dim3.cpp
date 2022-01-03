# include <cstdlib>
# include <sstream>
# include <string>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <mpi.h>

int main( int nargs, char* argv[] )
{
	MPI_Init(&nargs, &argv);
    int rank, nbp;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nbp);
    //~ std::cout << "Hello world from " 
              //~ << numero_du_processus << " in "
              //~ << nombre_de_processus << " executed" 
              //~ << std::endl;
    
    int tag = 78; // on aime bien avoir même tag dans tous les processus
              
    if (rank == 0)
    {
		MPI_Status status;
		int buf = 42;
		int last_rank = nbp-1; // car rank commence à 0
		MPI_Send(&buf, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
		//~ MPI_Recv(&buf, 1, MPI_INT, last_rank, tag, MPI_COMM_WORLD, &status);
		
		std::cout << "La valeur du jeton est : " << buf << " et je suis le processus numéro " << rank << std::endl;
	}
	
	else
	{
		int buf;
		MPI_Status status;
		MPI_Recv(&buf, 1, MPI_INT, (rank-1)%nbp, tag, MPI_COMM_WORLD, &status);
		std::cout << "La valeur du jeton est : " << buf << " et je suis le processus numéro " << rank << std::endl;

		//~ MPI_Send(&buf, 1, MPI_INT, (rank+1)%nbp, tag, MPI_COMM_WORLD);
		
		//~ std::cout << "Mon jeton vaut : " << buf << " et je suis le processus numéro " << rank << std::endl;
	}
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}


// Pistes : faire boucle sur la dimension de 0 jusqu'à la dimension actuelle
// chaque processus envoie à la dimension rank (la sienne) + i où i varie de 0 à la dimension actuelle
