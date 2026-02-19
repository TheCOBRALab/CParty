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
void prune_plist(std::vector<elem_prob_s> &p, std::vector<pf_t> &pu, std::vector<elem_prob_s> &pp,std::vector<elem_prob_s> &pl,std::vector<elem_prob_s> &plin, sparse_tree &tree,double gamma){

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
            pp.push_back(e);
            if(tree.tree[e.i].pair == e.j){
                plin.push_back(e);
            } else {
                pl.push_back(e);
            }  
        }
    
    }

}
void remake_structure(std::string &structure, sparse_tree &tree){
    cand_pos_t n = structure.length();
    std::vector<int> pairs;
    for (cand_pos_t j = 0; j < n; ++j) {
        if (structure[j] == '(') pairs.push_back(j);
        
        if (structure[j] == ')') {
            int i = pairs.back();
            pairs.pop_back();
            if(!tree.weakly_closed(i,j)){
                structure[i] = '[';
                structure[j] = ']';
            }
        }
    }
}

/**
 * @brief Register a candidate
 * @param i start
 * @param j end
 * @param e energy of candidate "M(i,j)"
 */
void register_candidate(std::vector<cand_list_t> &CL, cand_pos_t const &i, cand_pos_t const &j, pf_t const &e) {
    CL[j].emplace_back(cand_entry_t(i, e));
}

pf_t W_final_pf::compute_MEA(sparse_tree &tree,double gamma){
    std::string structure = std::string(n, '.');
    
    std::vector <elem_prob_s> p; // All elements with probabilities for pairs > cutoff
    std::vector <elem_prob_s> pp; // probability paired list
    std::vector <elem_prob_s> pl; // probability paired list excluding input structure
    std::vector <elem_prob_s> plin; // probability paired list of only input structure
    std::vector<pf_t> pu(n+1,1.0); // Probabilitiy unpaired list
    
    // Fill p with all pairs/probs > cutoff
    plist_from_probs(p,samples,num_samples,n,1e-4 / (1 + gamma));

    // Prune list to obtain unpaired probabilities as well as
    // to have pp store only pairs whose p(i,j) > pu[i] + pu[j] 
    prune_plist(p,pu,pp,pl,plin,tree,gamma);

    std::vector<pf_t> Mi(n+1);
    std::vector<pf_t> Mi1(n+1);
    std::vector<cand_list_t> CL(n+1);

    pf_t MEA = 0.;
    cand_pos_t index=0;
    for (cand_pos_t i = n; i > 0; --i) {
        Mi[i] = pu[i];
        if ((pl[index].i == i) && (pl[index].i > pl[index].j)) ++index;

        for (cand_pos_t j = i + 1; j <= n; ++j) {
            Mi[j] = Mi[j - 1] + pu[j];

            for (auto it = CL[j].begin(); CL[j].end() != it && it->first >= i; ++it) {
                cand_pos_t k = it->first;
                pf_t EA = it->second + Mi[k-1];
                Mi[j] = std::max(Mi[j], EA);
            }

            if ((pl[index].i == i) && (pl[index].j == j)) {
                pf_t EA = Mi1[j-1];

                EA += 2 * gamma * pl[index].p;

                if (Mi[j] < EA) {
                Mi[j] = EA;
                register_candidate(CL,i,j,EA);
                }
                ++index;
            }
        }
        Mi.swap(Mi1);
    }
    MEA = std::max(MEA, Mi1[n]);

    MEAdat bdat(structure,gamma,pl,pu,CL,Mi1);
    mea_backtrack(bdat, 1, n, 0);
    structure = bdat.structure;
    remake_structure(structure,tree);

    for(cand_pos_t i =0;i<(cand_pos_t) plin.size();++i){
        MEA+=2 * gamma * plin[i].p;
        structure[plin[i].i-1] = '(';
        structure[plin[i].j-1] = ')';
    }
    this->MEA_structure = structure;
    return MEA;
}

void mea_backtrack(MEAdat &bdat,cand_pos_t i,cand_pos_t j, int pair){

    int fail = 1;
    std::vector<pf_t> &Mi = bdat.Mi;
    std::vector<pf_t> &pu = bdat.pu;
    std::vector<cand_list_t> &CL = bdat.CL;

    if (pair) {
        bdat.structure[i - 1]  = '(';
        bdat.structure[j - 1]  = ')';
        i++;
        j--;
         /* We've done this before in MEA() but didn't keep the results */
        Mi[i-1] = 0;
        Mi[i] = pu[i];
        for (cand_pos_t k = i + 1; k <= j; k++) {
            Mi[k] = Mi[k - 1] + pu[k];
            for (auto it = CL[k].begin(); CL[k].end() != it && it->first >= i; ++it){
                pf_t EA = it->second + Mi[it->first - 1];
                Mi[k] = std::max(Mi[k], EA);
            }
        }
    }

    pf_t prec = std::numeric_limits<double>::epsilon() * Mi[j];
    /* Mi values are filled, do the backtrace */
    while ((j > i) && (Mi[j] <= (Mi[j - 1] + pu[j] + prec))) {
        bdat.structure[j - 1] = '.';
        --j;
    }

    for (auto it = CL[j].begin(); CL[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (Mi[j] <= it->second + Mi[k - 1] + prec) {
            if (k > i + 3) mea_backtrack(bdat, i, k-1, 0);

            mea_backtrack(bdat, k, j, 1);
            fail = 0;
        }
    }
    if (fail && j > i){
        std::cerr << "ERROR: backtrack failed for MEA()" << std::endl;
        exit(1);
    }
}

// If I find the PK-free MEA where I exclude base pairs that are part of the input. And then just add those base pairs in after. Doesn't that find the MEA?
