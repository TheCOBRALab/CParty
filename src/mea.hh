#ifndef MEA_HEADER
#define MEA_HEADER

#include <base_types.hh>
#include "part_func.hh"

#include <vector>

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

// void mea_backtrack(MEAdat &bdat,sparse_tree &tree,cand_pos_t i, cand_pos_t j, int pair);


#endif