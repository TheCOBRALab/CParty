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
/**
 * If I have an unpaired matrix, I can try to walk through one and get all pairings.
 * The second pass, I will do the same stack for pushing, but I will also have another for the fatgraph
 * If either stack is empty and you get to and opening base (op), you push it into the string. If the distance between the two bases of the same
 * type is filled with only unpaired relative to its parent, you skip adding it to the fatgraph. Popping from the stack does not affect this calculation
 * because as long as the stack has something inside that an internal could form from, popping won't affect it.
 * So when adding another op to the stack when the stack is not empty. I get the cp from the ptable and look at whether there is anything between i and ip and jp and j
 * This should ensure that I am always
*/
void generate_pt(std::string &structure, std::vector<int> &fres, std::vector<int> &up, int n){
   std::vector<int> paren;
   std::vector<int> sb;
   int count = 0;
   for(int j = 0;j<n;++j){
      if(structure[j] == '('){
          paren.push_back(j);
          count = 0;
      }
      else if(structure[j] == '['){
          sb.push_back(j);
          count = 0;
      }
      else if(structure[j] == ')'){
          int i = paren.back();
          fres[i] = j;
          fres[j] = i;
          paren.pop_back();
          count = 0;
      }
      else if(structure[j] == ']'){
          int i = sb.back();
          fres[i] = j;
          fres[j] = i;
          sb.pop_back();
          count = 0;
      }
      else{
        ++count;
        up[j] = count;
      }
   }
   if(!paren.empty() && !sb.empty()){
       std::cout << "Error: stacks aren't empty" << std::endl;
       exit(0);
   }
}
// i and j are the structured part
bool empty_region(std::vector<int> &up, int i, int j){
    if(up[j-1]>=j-i-1) return true;
    return false;
}

void generate_fatgraph(std::string &structure, std::string &fatgraph,std::vector<int> &fres, std::vector<int> &up,int n){
    std::vector<int> paren;
    std::vector<int> sb;
    std::string fatgraph_full = std::string(n,'.');
   for(int j = 0;j<n;++j){
      if(structure[j] == '('){
          if(paren.empty()){
            paren.push_back(j);
            fatgraph_full[j] ='(';
            fatgraph+='(';
          } else {
              int pparent = paren.back();
              if (!empty_region(up,pparent,j) || !empty_region(up,fres[j],fres[pparent])){
                paren.push_back(j);
                fatgraph_full[j] ='(';
                fatgraph+='(';
              }else{
                paren.push_back(j);
              }
          }
      }
      if(structure[j] == ')'){
          int i = paren.back();
          paren.pop_back();
          if(fatgraph_full[i] == '('){
              fatgraph_full[j] = ')';
              fatgraph+=')';
          }
      }

      if(structure[j] == '['){
          if(sb.empty()){
            sb.push_back(j);
            fatgraph_full[j] ='[';
            fatgraph+='[';
          } else {
              int pparent = sb.back();
              if (!empty_region(up,pparent,j) || !empty_region(up,fres[j],fres[pparent])){
                sb.push_back(j);
                fatgraph_full[j] ='[';
                fatgraph+='[';
              }else{
                sb.push_back(j);
              }
          }
      }
      if(structure[j] == ']'){
          int i = sb.back();
          sb.pop_back();
          if(fatgraph_full[i] == '['){
              fatgraph_full[j] = ']';
              fatgraph+=']';
          }
      }
   }
}

std::string W_final_pf::get_fatgraph(std::string structure){
    int n = structure.length();
    std::vector<int> fres;
    std::vector<int> up;
    fres.resize(n,-2);
    up.resize(n,0);
    generate_pt(structure,fres,up,n);
    std::string fatgraph = "";
    generate_fatgraph(structure,fatgraph,fres,up,n);
    return fatgraph;
}