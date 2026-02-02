#include "base_types.hh"
#include "part_func.hh"

#include <string>
#include <iostream>
#include <vector>


/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
std::string W_final_pf::compute_centroid(sparse_tree &tree, pf_t &dist, pf_t &diversity){
    dist = 0;
    diversity = 0;
    pf_t p = 0;
    std::string centroid = std::string(n, '.');

    for (cand_pos_t i = 1; i <= n; i++){
        for (cand_pos_t j = i + 1; j <= n; j++) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            p = (pf_t)samples[base_pair] / num_samples;
            diversity += p*(1.0-p);
            if (p > 0.5) {
                /* regular base pair */
                if(tree.weakly_closed(i,j)){
                    centroid[i - 1] = '(';
                    centroid[j - 1] = ')';
                }
                else{
                    centroid[i - 1] = '[';
                    centroid[j - 1] = ']';
                }
                dist += (1 - p);
            } else {
                dist += p;
            }
        }
    }
    diversity*=2; // As there are two sides of a base pair
    return centroid;
}


/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
std::string W_final_pf::compute_centroid_PK_only(sparse_tree &tree, pf_t &dist, pf_t &diversity){
    dist = 0;
    diversity = 0;
    pf_t p = 0;
    std::string centroid = std::string(n, '.');

    // I have to recalculate sample frequencies using only structures that contain pseudoknots
    std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples_PK;
    int num_samples_PK = 0;
    for(auto &it: structures){
        std::string structure = it.first;
        int frequency = it.second;
        if((structure.find('[') != std::string::npos)){
            cand_pos_t length = structure.length();
            num_samples_PK += frequency;
            std::vector<int> paren;
            std::vector<int> sb;
            for(cand_pos_t j=0;j<length;++j){
                if(structure[j] == '(') {
                    paren.push_back(j);
                    continue;
                }
                if(structure[j] == '[') {
                    sb.push_back(j);
                    continue;
                }

                if (structure[j] == ')'){
                    int x = paren[paren.size()-1];
                    paren.pop_back();
                    std::pair<cand_pos_tu, cand_pos_tu> base_pair(x+1,j+1);
                    samples_PK[base_pair]+=frequency;
                }
                if (structure[j] == ']'){
                    int x = sb[sb.size()-1];
                    sb.pop_back();
                    std::pair<cand_pos_tu, cand_pos_tu> base_pair(x+1,j+1);
                    samples_PK[base_pair]+=frequency;
                }
            }
        }
    }

    //Calculate centroid based on PK samples
    for (cand_pos_t i = 1; i <= n; i++){
        for (cand_pos_t j = i + 1; j <= n; j++) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            p = (pf_t)samples_PK[base_pair] / num_samples_PK;
            diversity += p*(1.0-p);
            if (p > 0.5) {
                /* regular base pair */
                if(tree.weakly_closed(i,j)){
                    centroid[i - 1] = '(';
                    centroid[j - 1] = ')';
                }
                else{
                    centroid[i - 1] = '[';
                    centroid[j - 1] = ']';
                }
                dist += (1 - p);
            } else {
                dist += p;
            }
        }
    }
    diversity*=2; // As there are two sides of a base pair
    return centroid;
}