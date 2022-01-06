# Mesure des temps (fichier fourier_mesure_temps)


Temps sélection des p% plus gros coefficients 0.0002542 secondes
Temps reconstitution de l'image compressée 0.0430296 secondes
Temps total compression image : 0.334036 secondes


# Parallélisation OMP (fichier fourier_compression_omp)

Pour la boucle :

for (std::uint32_t i = 1; i <= height; ++i )
        sparse.begin_rows[i] = sparse.begin_rows[i-1] + nbCoefsPerRow[i-1];

Il n'est pas possible de la paralléliser puisque pour créer l'élément i du tableau begin_rows, nous avons besoin
de l'élément précédent i-1, la construction doit donc se faire dans l'ordre

Ce sont donc les boucles for liées à la compression et la décompression qui peuvent être paralléliser avec OpenMP.
C'est une parallélisation sur les pixels qui est ici effectuée.

  Séquentiel      | time    | Speedup | 
------------------|---------|---------|
1                 |         |         | 
2                 |         |         |
3                 |         |         |
4                 |         |         | 
5                 |         |         |
6                 |         |         |
7                 |         |         |
8                 |         |         | 

En moyenne on trouve les temps suivants : 
Temps sélection des p% plus gros coefficients 0.0002847 secondes
Temps reconstitution de l'image compressée 0.0081761 secondes
Temps total compression image : 0.0897866 secondes

On a donc une accélération sur les différentes parties de :
Accélération sélection des p% plus gros coefficients : 0,89 -> NORMAL pas de parallélisation sur cette partie (donc environ = 1)
Accélération reconstitution de l'image compressée 5,26
Accélération totale compression image : 3,72
On a donc un temps de compression plus de 3 fois plus rapide grâce à la parallélisation OpenMP du fait des parallélisation des boucles for sur les pixels.

# Première parallélisation avec MPI (fichier fourier_compression_mpi_1_bis)

Problème sur la fonction Gather
Je ne trouve malheureusement pas l'erreur ce qui m'a fait perdre du temps et je ne peux pas aller plus loin
Je pense néanmoins que le code est quasiment bon, un détail m'échappe


