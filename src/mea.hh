#ifndef MEA_HEADER
#define MEA_HEADER

#include "base_types.hh"
#include "part_func.hh"

#include <vector>

typedef std::pair<cand_pos_t, pf_t> cand_entry_t;
typedef std::vector<cand_entry_t> cand_list_t;

/**
 *  @brief  Data structure which holds all relevant information for MEA probability entries
 */
struct elem_prob_s{
  cand_pos_t i;      /**<  @brief  Start position (usually 5' nucleotide that starts the element, e.g. base pair) */
  cand_pos_t j;      /**<  @brief  End position (usually 3' nucleotide that ends the element, e.g. base pair) */
  pf_t p;            /**<  @brief  Probability of the element */

    elem_prob_s(cand_pos_t start, cand_pos_t end, pf_t prob) {
        i = start;
        j = end;
        p = prob;
    }
};

struct MEAdat {
    std::vector<cand_pos_t> index;
    std::vector<elem_prob_s> pp;
    std::vector<elem_prob_s> plpk;
    std::vector<pf_t> pu;
    std::vector<pf_t> M;
    std::vector<pf_t> BE;
    std::vector<pf_t> WMBP;
    double gamma;
    std::vector<cand_list_t> CL;
    std::vector<cand_list_t> CLPK;
    std::string structure;
    MEAdat(std::vector<cand_pos_t> &index,std::vector<elem_prob_s> &pp,std::vector<elem_prob_s> &plpk,std::vector<pf_t> &pu,std::vector<pf_t> &M,std::vector<pf_t> &BE,std::vector<pf_t> &WMBP,double &gamma,std::vector<cand_list_t> &CL,std::vector<cand_list_t> &CLPK, std::string &structure){
        this->index = index;
        this->pp = pp;
        this->plpk = plpk;
        this->pu = pu;
        this->M = M;
        this->BE = BE;
        this->WMBP = WMBP;
        this->gamma = gamma;
        this->CL = CL;
        this->CLPK = CLPK;
        this->structure = structure;

    }
};

struct Cand_comp {
        bool operator()(const cand_entry_t &x, cand_pos_t y) const { return x.first > y; }
} cand_comp;

/*
 *  sort by sequence position:
 *  1. in descending order for i
 *  2. in ascending order for j
 */
int comp_plist(const elem_prob_s &a,const elem_prob_s &b){
    // cand_pos_t di  = (b.i - a.i);
    // if (di != 0) return di;
    // return a.j - b.j;
    if(a.i != b.i) return a.i > b.i;
    return a.j < b.j;
}

void plist_from_probs(std::vector<elem_prob_s> &p,std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples, cand_pos_t n, int num_samples, double cutoff);

void prune_plist(std::vector<elem_prob_s> &p, std::vector<pf_t> &pu, std::vector<elem_prob_s> &pl, double gamma);

void mea_backtrack_pk(MEAdat &bdat,sparse_tree &tree,cand_pos_t i,cand_pos_t j, pf_t e);
void mea_backtrack(MEAdat &bdat,sparse_tree &tree,cand_pos_t i,cand_pos_t j, int pair);


#endif