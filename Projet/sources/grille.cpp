#include "grille.hpp"

épidémie::Grille::Grille( std::size_t population )
{
    int taille = int(sqrt(population));
    m_dim_x = (population+taille-1)/taille;
    m_dim_y = taille;
    std::vector<StatistiqueParCase>(m_dim_x*m_dim_y).swap(m_statistiques);
}
//
std::size_t 
épidémie::Grille::nombreTotalContaminésGrippe() const
{
    std::size_t total = 0;
    for (auto statistique : m_statistiques)
    {
        total += statistique.nombre_contaminant_grippé_et_contaminé_par_agent;
        total += statistique.nombre_contaminant_seulement_grippé;
    }
    return total;
}
//
std::size_t 
épidémie::Grille::nombreTotalContaminésAgentPathogène() const
{
    std::size_t total = 0;
    for (auto statistique : m_statistiques )
    {
        total += statistique.nombre_contaminant_seulement_contaminé_par_agent;
        total += statistique.nombre_contaminant_grippé_et_contaminé_par_agent;
    }
    return total;
}
//

std::vector<int>
épidémie::Grille::statistiquesToVector() const
{
    // Allocate vec for the three ints that are in each statistic
    std::vector<int> vec(3*m_statistiques.size());
    for (auto i = 0; i < m_statistiques.size(); i++)
    {
        vec[i*3] = m_statistiques[i].nombre_contaminant_seulement_grippé;
        vec[i*3 + 1] = m_statistiques[i].nombre_contaminant_seulement_contaminé_par_agent;
        vec[i*3 + 2] = m_statistiques[i].nombre_contaminant_grippé_et_contaminé_par_agent;
    }
    return vec;
}
//
void
épidémie::Grille::vectorToStatistiques(std::vector<int> vec)
{
    size_t num_statistiques = vec.size()/3;
    if (num_statistiques != m_statistiques.size())
        m_statistiques.resize(num_statistiques);
    
    for (auto i = 0; i < num_statistiques; i++) {
        m_statistiques[i].nombre_contaminant_seulement_grippé = vec[i*3];
        m_statistiques[i].nombre_contaminant_seulement_contaminé_par_agent = vec[i*3 + 1];
        m_statistiques[i].nombre_contaminant_grippé_et_contaminé_par_agent = vec[i*3 + 2];
    }
}
