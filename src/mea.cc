#include "mea.hh"
#include "mwm.hh"

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#define debug 0

int compute_exterior_cases(cand_pos_t l, cand_pos_t j, sparse_tree &tree) {
    // Case 1 -> l is not covered
    bool case1 = tree.tree[l].parent->index <= 0;
    // Case 2 -> l is paired
    bool case2 = tree.tree[l].pair > 0;
    // Case 3 -> l is part of a closed subregion
    // bool case3 = 0;
    // Case 4 -> l.bp(l) i.e. l.j does not cross anything -- could I compare parents instead?
    bool case4 = j < tree.Bp(l, j);
    // By bitshifting each one, we have a more granular idea of what cases fail and is faster than branching
    return (case1 << 2) | (case2 << 1) | case4;
}

/**
 * 
 * @brief Given the probabilities found prior, fill the vector p of all entries whose value is greater than the cutoff
 * 
 */
void plist_from_probs(std::vector<elem_prob_s> &p,std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples, int num_samples, cand_pos_t n, double cutoff){
    for (cand_pos_t i = n; i >= 1; --i) {
        for (cand_pos_t j = i; j <= n; ++j) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            pf_t prob = (pf_t)samples[base_pair] / num_samples;
            if(prob < cutoff) continue;

            elem_prob_s s(i,j,prob);
            p.push_back(s);

        }
    }
}

  /*
   * produce a list containing all base pairs with
   * 2*gamma*p_ij > p^u_i + p^u_j
   * already sorted to be in the order we need them within the DP
   */
void prune_plist(std::vector<elem_prob_s> &p, std::vector<pf_t> &pu, std::vector<elem_prob_s> &pl, std::vector<elem_prob_s> &pp,  std::vector<elem_prob_s> &plpk, sparse_tree& tree,double gamma){

    unsigned int nt = 2; // Two sides of base pair
    pf_t pug; // prob unpaired gamma

    // For me, it's always the case that I am looking at a base pair. There is no chance of gquad or anything
    // Subtracting base pair to find probability of unpaired
    for(elem_prob_s pc: p){
        pu[pc.i] -= pc.p;
        pu[pc.j] -= pc.p;
    }

    // This is important because I am saying here that I already know that all i.j I have are better than i + j
    for(elem_prob_s e: p){

        pug = pu[e.i] + pu[e.j];

        if (e.p * nt * gamma > pug) {
            pl.push_back(e);
            if(tree.weakly_closed(e.i,e.j)){
                pp.push_back(e);
            } else {
                plpk.push_back(e);
            }  
        }
    
    }

}

pf_t W_final_pf::compute_MEA(sparse_tree &tree,double gamma){
    std::string structure = std::string(n, '.');
    
    std::vector <elem_prob_s> p; // All elements with probabilities for pairs > cutoff
    std::vector <elem_prob_s> plpk; // probability paired list PK
    std::vector <elem_prob_s> pp; // probability paired list
    std::vector <elem_prob_s> pl;
    std::vector<pf_t> pu(n+1,1.0); // Probabilitiy unpaired list
    
    // Fill p with all pairs/probs > cutoff
    plist_from_probs(p,samples,num_samples,n,1e-4 / (1 + gamma));

    // // Prune list to only those ...
    prune_plist(p,pu,pl,pp,plpk,tree,gamma);
    
    pf_t MEA = 0.0;

    maxWeightMatching MWW(pl,n+1);
    MEA = MWW.maximumWeightedMatching(pu,structure,tree);
    
    this->MEA_structure = structure;

    return MEA;
}