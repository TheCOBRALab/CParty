#include "mea.hh"

#include <string>
#include <iostream>
#include <vector>

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
void prune_plist(std::vector<elem_prob_s> &p, std::vector<pf_t> &pu, std::vector<elem_prob_s> &pl, double gamma){

    unsigned int nt;
    pf_t pug; // prob unpaired gamma

    // For me, it's always the case that I am looking at a base pair. There is no chance of gquad or anything
    // Subtracting base pair to find probability of unpaired
    for(elem_prob_s pc: p){
        pu[pc.i] -= pc.p;
        pu[pc.j] -= pc.p;
    }

    // This is important because I am saying here that I already know that all i.j I have are better than i + j
    for(elem_prob_s pc: p){

        nt  = 2;
        pug = pu[pc.i] + pu[pc.j];

        if (pc.p * nt * gamma > pug) {  
            pl.push_back(pc);
        }
    
    }

}

struct BlossomMatching {
    cand_pos_t n;
    std::vector<std::vector<cand_pos_t>> CL;
    std::vector<elem_prob_s> pl;

    std::vector<cand_pos_t> mate, parent, base;
    std::vector<cand_pos_t> type;
    std::vector<bool> inq, inb;
    std::vector<cand_pos_t> q;

    BlossomMatching(cand_pos_t n, std::vector<elem_prob_s> &pl) : n(n), CL(n), mate(n,-1), parent(n,-1), base(n), type(n), inq(n), inb(n) {
        this->pl = pl;
    }

    int lca(int a, int b) {
        static std::vector<bool> used;
        used.resize(n, false);
        while (true) {
            a = base[a];
            used[a] = true;
            if (mate[a] == -1) break;
            a = parent[mate[a]];
        }
        while (true) {
            b = base[b];
            if (used[b]) return b;
            b = parent[mate[b]];
        }
    }

    void contract_blossom(cand_pos_t v, cand_pos_t u, cand_pos_t b, std::queue<cand_pos_t> Q){
        auto mark = [&](int x, int lca) {
            while (base[x] != lca) {
                int m = mate[x];
                int p = parent[m];
                base[x] = base[m] = lca;
                if (type[m] == 1) { // inner
                    type[m] = 0;    // make outer
                    Q.push(m);
                }
                x = p;
            }
        };

        mark(v, b);
        mark(u, b);
    }

    bool bfs(cand_pos_t root) {
        fill(parent.begin(), parent.end(), -1);
        fill(type.begin(), type.end(), -1);
        for(cand_pos_t i = 0;i<(cand_pos_t) base.size();++i) base[i] = i;

        std::queue<cand_pos_t> Q;
        Q.push(root);
        type[root] = 0;        // EVEN (outer)

        while (!Q.empty()) {
            cand_pos_t v = Q.front(); Q.pop();

            for (cand_pos_t eid : CL[v]) {
                const elem_prob_s &e = pl[eid];
                cand_pos_t u = (e.i == v ? e.j : e.i);

                // Skip invalid transitions
                if (base[u] == base[v] || mate[v] == u)
                    continue;

                if (type[u] == -1) {
                    // Unvisited vertex
                    type[u] = 1;      // ODD (inner)
                    parent[u] = v;

                    if (mate[u] == -1) {
                        // Augmenting path found!
                        augment(u);
                        return true;
                    }

                    cand_pos_t m = mate[u];
                    type[m] = 0;      // EVEN (outer)
                    Q.push(m);

                } else if (type[u] == 0) {
                    // Found EVEN → EVEN: possible blossom
                    cand_pos_t b = lca(v, u);
                    contract_blossom(v, u, b, Q);
                }
            }
        }
        return false;
    }

    void augment(cand_pos_t u){
        while(u!=-1){
            cand_pos_t p = parent[u];
            cand_pos_t pi = (p == -1 ? -1 : mate[p]);

            mate[u]=p;

            if(p!=-1) mate[p]=u;
            
            u=pi;
        }
    }

    pf_t max_weight_matching(std::vector<pf_t> &pu,std::string &structure, sparse_tree &tree){
        // Greedy
        std::sort(pl.begin(),pl.end(),[](const elem_prob_s &a, const elem_prob_s &b) {return a.p > b.p;});

        for(auto &e: pl){
            if(mate[e.i]==-1 && mate[e.j]==-1){
                mate[e.i]=e.j;
                mate[e.j]=e.i;
            }
        }

        for(cand_pos_t i=0;i<n;i++){
            if(mate[i]==-1){
                cand_pos_t u = bfs(i);
                if(u!=-1) augment(u);
            }
        }

        pf_t MEA=0.0;
        for(elem_prob_s &e: pl){
            if(mate[e.i]==e.j) {
                MEA+=2*e.p;
                if(tree.weakly_closed(e.i,e.j)){
                    structure[e.i-1] = '(';
                    structure[e.j-1] = ')';
                } else{
                    structure[e.i-1] = '[';
                    structure[e.j-1] = ']';
                }
            }
        }
        for(cand_pos_t i = n; i>0;--i) if(mate[i] == -1) MEA += pu[i];
        return MEA;
    }

};

pf_t W_final_pf::compute_MEA(sparse_tree &tree,double gamma){
    std::string structure = std::string(n, '.');
    
    std::vector <elem_prob_s> p; // All elements with probabilities for pairs > cutoff
    std::vector <elem_prob_s> pl; // probability list
    std::vector <elem_prob_s> pp; // probability paired list
    std::vector<pf_t> pu; // Probabilitiy unpaired list
    pu.resize(n+1,1.0);
    
    // Fill p with all pairs/probs > cutoff
    plist_from_probs(p,samples,num_samples,n,1e-4 / (1 + gamma));

    // // Prune list to only those ...
    prune_plist(p,pu,pl,gamma);

    // // sort by i then by j and copy this to pp
    std::sort(pl.begin(),pl.end(),comp_plist);
    pp = pl;
    
    pf_t MEA = 0.0;

    BlossomMatching BM(n+1,pl);
    MEA = BM.max_weight_matching(pu,structure,tree);
    this->MEA_structure = structure;

    return MEA;
}