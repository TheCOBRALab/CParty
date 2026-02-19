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

void plist_from_probs(std::vector<elem_prob_s> &p,std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples, cand_pos_t n, int num_samples, double cutoff);

void prune_plist(std::vector<elem_prob_s> &p, std::vector<pf_t> &pu, std::vector<elem_prob_s> &pl, std::vector<elem_prob_s> &pp,  std::vector<elem_prob_s> &plpk, sparse_tree& tree,double gamma);


#endif