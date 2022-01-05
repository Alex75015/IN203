#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include "contexte.hpp"
#include "individu.hpp"
#include "graphisme/src/SDL2/sdl2.hpp"
#include <chrono>
# include <mpi.h>


MPI_Comm simulationComm;

// graine var globale

void màjStatistique( épidémie::Grille& grille, std::vector<épidémie::Individu> const& individus )
{
    for ( auto& statistique : grille.getStatistiques() )
    {
        statistique.nombre_contaminant_grippé_et_contaminé_par_agent = 0;
        statistique.nombre_contaminant_seulement_contaminé_par_agent = 0;
        statistique.nombre_contaminant_seulement_grippé              = 0;
    }
    auto [largeur,hauteur] = grille.dimension();
    auto& statistiques = grille.getStatistiques();
    for ( auto const& personne : individus )
    {
        auto pos = personne.position();

        std::size_t index = pos.x + pos.y * largeur;
        if (personne.aGrippeContagieuse() )
        {
            if (personne.aAgentPathogèneContagieux())
            {
                statistiques[index].nombre_contaminant_grippé_et_contaminé_par_agent += 1;
            }
            else 
            {
                statistiques[index].nombre_contaminant_seulement_grippé += 1;
            }
        }
        else
        {
            if (personne.aAgentPathogèneContagieux())
            {
                statistiques[index].nombre_contaminant_seulement_contaminé_par_agent += 1;
            }
        }
    }
}

void afficheSimulation(sdl2::window& écran, épidémie::Grille const& grille, std::size_t jour)
{
    auto [largeur_écran,hauteur_écran] = écran.dimensions();
    auto [largeur_grille,hauteur_grille] = grille.dimension();
    auto const& statistiques = grille.getStatistiques();
    sdl2::font fonte_texte("./graphisme/src/data/Lato-Thin.ttf", 18);
    écran.cls({0x00,0x00,0x00});
    // Affichage de la grille :
    std::uint16_t stepX = largeur_écran/largeur_grille;
    unsigned short stepY = (hauteur_écran-50)/hauteur_grille;
    double factor = 255./15.;

    for ( unsigned short i = 0; i < largeur_grille; ++i )
    {
        for (unsigned short j = 0; j < hauteur_grille; ++j )
        {
            auto const& stat = statistiques[i+j*largeur_grille];
            int valueGrippe = stat.nombre_contaminant_grippé_et_contaminé_par_agent+stat.nombre_contaminant_seulement_grippé;
            int valueAgent  = stat.nombre_contaminant_grippé_et_contaminé_par_agent+stat.nombre_contaminant_seulement_contaminé_par_agent;
            std::uint16_t origx = i*stepX;
            std::uint16_t origy = j*stepY;
            std::uint8_t red = valueGrippe > 0 ? 127+std::uint8_t(std::min(128., 0.5*factor*valueGrippe)) : 0;
            std::uint8_t green = std::uint8_t(std::min(255., factor*valueAgent));
            std::uint8_t blue= std::uint8_t(std::min(255., factor*valueAgent ));
            écran << sdl2::rectangle({origx,origy}, {stepX,stepY}, {red, green,blue}, true);
        }
    }

    écran << sdl2::texte("Carte population grippée", fonte_texte, écran, {0xFF,0xFF,0xFF,0xFF}).at(largeur_écran/2, hauteur_écran-20);
    écran << sdl2::texte(std::string("Jour : ") + std::to_string(jour), fonte_texte, écran, {0xFF,0xFF,0xFF,0xFF}).at(0,hauteur_écran-20);
    écran << sdl2::flush;
}

void Processus_Affichage(bool affiche, épidémie::ContexteGlobal contexte)
{
	constexpr const unsigned int largeur_écran = 1280, hauteur_écran = 1024;
	sdl2::window écran("Simulation épidémie de grippe", {largeur_écran,hauteur_écran});
	
    // contexte.déplacement_maximal = 1; <= Si on veut moins de brassage
    // contexte.taux_population = 400'000;
    //contexte.taux_population = 1'000;

    épidémie::Grille grille{contexte.taux_population};
    
    std::vector<int> statVector = grille.statistiquesToVector();
	
	
	sdl2::event_queue queue;
	bool quitting = false;
	
	MPI_Request terminer_requete;
	int jours_écoulés;
	
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //~ std::chrono::time_point < std::chrono::system_clock > start, end;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//~ start = std::chrono::system_clock::now();	
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	
	while (true)
    {
		auto events = queue.pull_events();
        for ( const auto& e : events)
        {
            if (e->kind_of_event() == sdl2::event::quit)
            {
                quitting = true;
                MPI_Isend( &quitting , 1 , MPI_INT , 1 , 0 , MPI_COMM_WORLD, &terminer_requete);
            }
        }
        
        if (quitting)
            break;
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//~ start = std::chrono::system_clock::now();	
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	
		
        //#############################################################################################################
        //##    Affichage des résultats pour le temps  actuel
        //#############################################################################################################
        
        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//~ start = std::chrono::system_clock::now();	
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
		MPI_Status status;
		
		if (affiche)
		{
			MPI_Request wait_request;
            
            MPI_Isend( nullptr , 0 , MPI_INT , 1 , 1 , MPI_COMM_WORLD , &wait_request);
            MPI_Recv( statVector.data() , statVector.size() , MPI_INT , 1 , MPI_ANY_TAG , MPI_COMM_WORLD , &status); 
            grille.vectorToStatistiques(statVector);
        
            jours_écoulés = status.MPI_TAG;
            afficheSimulation(écran, grille, jours_écoulés);
        }   
		
	}// Fin boucle temporelle
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	int flag = 0;
	MPI_Iprobe( 1 , 0 , MPI_COMM_WORLD , &flag , MPI_STATUS_IGNORE);
		
	if (!flag && affiche) {
	MPI_Recv( statVector.data() , statVector.size() , MPI_INT , 1 , MPI_ANY_TAG , 
				MPI_COMM_WORLD , MPI_STATUS_IGNORE);
	}
        
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//~ end = std::chrono::system_clock::now();
	//~ std::chrono::duration < double >elapsed_seconds = end - start;
	//~ std::cout << "Temps  : " << elapsed_seconds.count() << " secondes\n";
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	/*std::cout << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
			  << grille.nombreTotalContaminésAgentPathogène() << std::endl;*/
	
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//~ end = std::chrono::system_clock::now();
	//~ std::chrono::duration < double >elapsed_seconds = end - start;
	//~ std::cout << "Temps  : " << elapsed_seconds.count() << " secondes\n";
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	
    
    //~ if (affiche)
	//~ {
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//~ end = std::chrono::system_clock::now();
	//~ std::chrono::duration < double >elapsed_seconds = end - start;
	//~ std::cout << "Temps moyen par pas de temps à l'affichage : " << elapsed_seconds.count()/jours_écoulés << " secondes\n";
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//~ }
}

void Processus_Simulation(bool affiche, épidémie::ContexteGlobal contexte)
{
    int nbp;
	MPI_Comm_size(simulationComm, &nbp);

	int rank;
	MPI_Comm_rank(simulationComm, &rank);
	
	
	int qtd_population = contexte.taux_population / nbp + (rank < contexte.taux_population % nbp);
    int i_begin = rank * contexte.taux_population / nbp + (rank < contexte.taux_population % nbp ? rank : contexte.taux_population % nbp);
    int i_end = i_begin + qtd_population;
    
	unsigned int graine_aléatoire = 1;
	std::uniform_real_distribution<double> porteur_pathogène(0.,1.);
    // contexte.déplacement_maximal = 1; <= Si on veut moins de brassage
    // contexte.taux_population = 400'000;
    //contexte.taux_population = 1'000;

    std::vector<épidémie::Individu> population;
    population.reserve(qtd_population);
    épidémie::Grille grille{contexte.taux_population};

    auto [largeur_grille,hauteur_grille] = grille.dimension();
    // L'agent pathogène n'évolue pas et reste donc constant...
    épidémie::AgentPathogène agent(graine_aléatoire++);
    // Initialisation de la population initiale :
    
    graine_aléatoire += i_begin;
    for (std::size_t i = i_begin; i < i_end; ++i )
    {
        std::default_random_engine motor(100*(i+1));
        population.emplace_back(graine_aléatoire++, contexte.espérance_de_vie, contexte.déplacement_maximal);
        population.back().setPosition(largeur_grille, hauteur_grille);
        if (porteur_pathogène(motor) < 0.2)
        {
            population.back().estContaminé(agent);   
        }
    }

    std::ofstream output;
    if (rank == 0) {
        output.open("Courbe.dat");
		output << "# jours_écoulés \t nombreTotalContaminésGrippe \t nombreTotalContaminésAgentPathogène()" << std::endl;
	}
	
	std::size_t jours_écoulés = 0;
    int         jour_apparition_grippe = 0;
    int         nombre_immunisés_grippe= (contexte.taux_population*23)/100;
    
	épidémie::Grippe grippe(0);
    //~ bool quitting = false; // fonctionne pas pour MPI_Probe donc utilise un int
    int quitting = 0;
	
    
    
    std::chrono::time_point < std::chrono::system_clock > begin, end;
    if (rank == 0){
		MPI_Iprobe( 0, 0 , MPI_COMM_WORLD , &quitting , MPI_STATUS_IGNORE);

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

		std::cout << "Début boucle épidémie" << std::endl << std::flush;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		begin = std::chrono::system_clock::now();
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	}
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    while (!quitting)
    {
		
		//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
		//~ MPI_Iprobe( 0, 0 , MPI_COMM_WORLD , &quitting , MPI_STATUS_IGNORE);
		//~ if (quitting)
            //~ break;


        if (jours_écoulés%365 == 0)// Si le premier Octobre (début de l'année pour l'épidémie ;-) )
        {
            grippe = épidémie::Grippe(jours_écoulés/365);
            jour_apparition_grippe = grippe.dateCalculImportationGrippe();
            grippe.calculNouveauTauxTransmission();
            // 23% des gens sont immunisés. On prend les 23% premiers
            for ( int ipersonne = i_begin; ipersonne < nombre_immunisés_grippe && ipersonne < i_end; ++ipersonne)
            {
                population[ipersonne - i_begin].devientImmuniséGrippe();
            }
            for ( int ipersonne = nombre_immunisés_grippe > i_begin ? nombre_immunisés_grippe : i_begin; ipersonne < i_end; ++ipersonne )
            {
                population[ipersonne - i_begin].redevientSensibleGrippe();
            }
        }
        if (jours_écoulés%365 == std::size_t(jour_apparition_grippe))
        {
            for (int ipersonne = nombre_immunisés_grippe > i_begin ? nombre_immunisés_grippe : i_begin;
                 ipersonne < nombre_immunisés_grippe + 25 && ipersonne < i_end; ++ipersonne )
            {
                population[ipersonne - i_begin].estContaminé(grippe);
            }
        }
        // Mise à jour des statistiques pour les cases de la grille :
        màjStatistique(grille, population);
        // On parcout la population pour voir qui est contaminé et qui ne l'est pas, d'abord pour la grippe puis pour l'agent pathogène
        
        std::vector<int> statVector = grille.statistiquesToVector();
        MPI_Allreduce( MPI_IN_PLACE , statVector.data(), statVector.size() , MPI_INT , MPI_SUM , simulationComm);
        grille.vectorToStatistiques(statVector);
        
        std::size_t compteur_grippe = 0, compteur_agent = 0, mouru = 0;
        for ( auto& personne : population )
        {
            if (personne.testContaminationGrippe(grille, contexte.interactions, grippe, agent))
            {
                compteur_grippe ++;
                personne.estContaminé(grippe);
            }
            if (personne.testContaminationAgent(grille, agent))
            {
                compteur_agent ++;
                personne.estContaminé(agent);
            }
            // On vérifie si il n'y a pas de personne qui veillissent de veillesse et on génère une nouvelle personne si c'est le cas.
            if (personne.doitMourir())
            {
                mouru++;
                unsigned nouvelle_graine = jours_écoulés + personne.position().x*personne.position().y;
                personne = épidémie::Individu(nouvelle_graine, contexte.espérance_de_vie, contexte.déplacement_maximal);
                personne.setPosition(largeur_grille, hauteur_grille);
            }
            personne.veillirDUnJour();
            personne.seDéplace(grille);
        }
        //~ MPI_Iprobe( 0, 0 , MPI_COMM_WORLD , &quitting , MPI_STATUS_IGNORE);
        //~ if (quitting)
            //~ break;
        
        
        if (rank == 0 && affiche)  
        {
            if (affiche) {
                int waiting;
                
                MPI_Iprobe( 0, 1 , MPI_COMM_WORLD , &waiting , MPI_STATUS_IGNORE);
                
                if (waiting) 
                {
                    std::vector<int> statVector = grille.statistiquesToVector();
                    MPI_Recv( nullptr , 0 , MPI_INT , 0 , 1 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
                    MPI_Send( statVector.data() , statVector.size() , MPI_INT , 0 , jours_écoulés , MPI_COMM_WORLD);
                }
            }
        }
        //~ MPI_Send(&jours_écoulés, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	


        /*std::cout << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
                  << grille.nombreTotalContaminésAgentPathogène() << std::endl;*/

	if (rank == 0) 
	{
		output << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
			   << grille.nombreTotalContaminésAgentPathogène() << std::endl;

		// On regarde si il faut quitter
		MPI_Iprobe( 0, 0 , MPI_COMM_WORLD , &quitting , MPI_STATUS_IGNORE);
	}

	MPI_Bcast( &quitting , 1 , MPI_INT , 0 , simulationComm);

	jours_écoulés += 1;
        
    }// Fin boucle temporelle
    
    if (rank == 0) 
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration < double >elapsed_seconds = end - begin;
        std::cout << "Temps moyen par pas de temps à l'affichage : "
                        << elapsed_seconds.count()/jours_écoulés 
                        << " secondes" << std::endl;
    
        
        MPI_Request terminer_requete;
        MPI_Isend( nullptr , 0 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , &terminer_requete);

        output.close();
    }
}

void simulation(bool affiche)
{
    
    int nbp;
	MPI_Comm_size(MPI_COMM_WORLD, &nbp);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	épidémie::ContexteGlobal contexte;
	contexte.interactions.β = 60.;
	
	MPI_Comm_split( MPI_COMM_WORLD , rank != 0 , rank == 0 ? 0 : rank - 1 , &simulationComm);
    
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	if (rank==0)
	{
		Processus_Affichage(affiche, contexte);
	}
	
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	else
	{
		Processus_Simulation(affiche, contexte);
	}
	std::cout << "Le processus" << rank << " est terminé\n";
}

int main(int argc, char* argv[])
{
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	
    // parse command-line
    bool affiche = true;
    for (int i=0; i<argc; i++) {
      std::cout << i << " " << argv[i] << "\n";
      if (std::string("-nw") == argv[i]) affiche = false;
    }
  
    sdl2::init(argc, argv);
    {
        MPI_Init( &argc, &argv );
        simulation(affiche);
        MPI_Finalize();
    }
    sdl2::finalize();
    return EXIT_SUCCESS;
}

