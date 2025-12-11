#include "mea.hh"

#include <string>
#include <iostream>
#include <vector>

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
void prune_plist(std::vector<elem_prob_s> &p, std::vector<pf_t> &pu, std::vector<elem_prob_s> &pp,  std::vector<elem_prob_s> &plpk, sparse_tree& tree,double gamma){

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
            if(tree.weakly_closed(e.i,e.j)){
                pp.push_back(e);
            } else {
                plpk.push_back(e);
            }  
        }
    
    }

}

struct BlossomMatching {
    cand_pos_t n;
    std::vector<std::vector<cand_pos_t>> CL;
    std::vector<elem_prob_s> pl;

    std::vector<cand_pos_t> mate, parent, base;
    std::vector<pf_t> path_weight;
    std::vector<cand_pos_t> type;
    std::vector<cand_pos_t> q;

    BlossomMatching(cand_pos_t n, std::vector<elem_prob_s> &pl) : n(n), CL(n), mate(n,-1), parent(n,-1), base(n), path_weight(n), type(n) {
        this->pl = pl;
        // We sort by descending probability for when we do a greedy addition of base pairs. This is also better for filling the candidate list
        std::sort(pl.begin(),pl.end(),[](const elem_prob_s &a, const elem_prob_s &b) {return a.p > b.p;});

        for (cand_pos_t eid = 0; eid < (cand_pos_t) pl.size(); ++eid) {
            const elem_prob_s &e = pl[eid];
            CL[e.i].push_back(eid);
            CL[e.j].push_back(eid);
        }

    }

    cand_pos_t lca(cand_pos_t a, cand_pos_t b) {
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

    void contract_blossom(cand_pos_t v, cand_pos_t u, cand_pos_t b, std::queue<cand_pos_t> &Q){
        auto mark = [&](cand_pos_t x, cand_pos_t lca) {
            while (base[x] != lca) {
                cand_pos_t m = mate[x];
                cand_pos_t p = parent[m];
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

    //
    void augment(cand_pos_t u){
        // I know the path from bfs for pairings
        // u will be -1 when I get to the root
        while(u!=-1){
            cand_pos_t p = parent[u]; // Get the base pair index for u i.e the parent.
            cand_pos_t pi = (p == -1 ? -1 : mate[p]); // If we are at the root, 

            mate[u]=p;

            if(p!=-1) mate[p]=u;
            
            u=pi;
        }
    }

    /** 
     * @brief Breadth first search from the current root node looking for augmentable paths
     * @param parent Defines the parents for each vertex
     * @param type Defines the type of path we are looking at. -1 is unvisited, 0 is even distance from root, 1 is odd distance from root
     * @param base
     * @param mate List which defines the current pairings
     * @param CL Candidate list or adjacency list which gives the edge list (pl) index for pairing for an index
     * 
    */
    void bfs(cand_pos_t root) {
        // Reset parent, type, and base to default values
        fill(parent.begin(), parent.end(), -1);
        fill(type.begin(), type.end(), -1);
        fill(path_weight.begin(), path_weight.end(), 0.0);
        for(cand_pos_t i = 0;i<(cand_pos_t) base.size();++i) base[i] = i;

        pf_t best_augment_weight      = 0.0;
        cand_pos_t best_augment_node  = -1;

        // Push root into queue and set it as an even path
        std::queue<cand_pos_t> Q;
        Q.push(root);
        type[root] = 0;

        while (!Q.empty()) {
            // Pop index from queue
            cand_pos_t v = Q.front(); Q.pop();

            // go through all of its possible pairings
            for (cand_pos_t eid : CL[v]) {
                // I pushed the edge index to both so I don't know if I'm looking at i or j in it
                const elem_prob_s &e = pl[eid];
                cand_pos_t u = (e.i == v ? e.j : e.i);

                // The first check ensures we are not inside a blossom. The second check is because we don't do anything for an already paired edge that we've looked at
                if (base[u] == base[v] || mate[v] == u) continue;

                pf_t candidate_weight = path_weight[v] + e.p;
                
                // We are at an unvisited node
                if (type[u] == -1 || candidate_weight> path_weight[u]) {
                    // This means that we are at the odd portion distance wise at the moment and v is the parent as u is one of the pairable indices to v.
                    type[u] = 1;
                    parent[u] = v;
                    path_weight[u] = candidate_weight;

                    // If u is unpaired currently within the path (odd -> even -> odd ...), we are at an augmentable path
                    if (mate[u] == -1) {
                        if(candidate_weight> best_augment_weight){
                            best_augment_node = u;
                            best_augment_weight = candidate_weight;
                        }
                        continue;
                    }

                    cand_pos_t m = mate[u];
                    path_weight[m] = candidate_weight;
                    type[m] = 0;
                    Q.push(m);

                } else if (type[u] == 0) {
                    // Found EVEN → EVEN: possible blossom
                    cand_pos_t b = lca(v, u);
                    contract_blossom(v, u, b, Q);
                }
            }
        }

        if (best_augment_node != -1)
        augment(best_augment_node);
    }

    pf_t max_weight_matching(std::vector<pf_t> &pu,std::string &structure, sparse_tree &tree){
        // We start by doing a greedy filling of the base pairs. This will reduce calls to bfs and augment later
        for(auto &e: pl){
            if(mate[e.i]==-1 && mate[e.j]==-1){
                mate[e.i]=e.j;
                mate[e.j]=e.i;
            }
        }

        for(cand_pos_t i=0;i<n;i++){
            // If it's a strictly unpaired base, skip
            if(!CL[i].empty()){
                // If we are not looking at an augmentable path, then skip
                if(mate[i]==-1){
                    bfs(i);
                }
            }
        }
        
        pf_t MEA=0.0;
        // As we now have the mate array/pairing array, we can calculate the structure and MEA value
        // simply by going through and changing the structure base pairs and adding the probabilities
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
        // If mate is -1, then it's an unpaired base and we add that probability
        for(cand_pos_t i = n; i>0;--i) if(mate[i] == -1) MEA += pu[i];
        return MEA;
    }

};

/**
 * @brief Register a candidate
 * @param i start
 * @param j end
 * @param e energy of candidate "V(i,j)"
 */
void register_candidate(std::vector<cand_list_t> &CL, cand_pos_t const &i, cand_pos_t const &j, pf_t const &e) {
    CL[j].emplace_back(cand_entry_t(i, e));
}

cand_pos_t get_index(std::vector<cand_pos_t> &index, cand_pos_t i, cand_pos_t j) {
        if(j<i) return 0; // index 0 will be 0
        return index[i] + j - i;
}

pf_t W_final_pf::compute_MEA(sparse_tree &tree,double gamma){
    std::string structure = std::string(n, '.');
    
    std::vector <elem_prob_s> p; // All elements with probabilities for pairs > cutoff
    std::vector <elem_prob_s> plpk; // probability paired list PK
    std::vector <elem_prob_s> pp; // probability paired list
    std::vector<pf_t> pu; // Probabilitiy unpaired list
    pu.resize(n+1,1.0);
    
    // Fill p with all pairs/probs > cutoff
    plist_from_probs(p,samples,num_samples,n,1e-4 / (1 + gamma));

    // // Prune list to only those ...
    prune_plist(p,pu,pp,plpk,tree,gamma);
    
    pf_t MEA = 0.0;

    std::vector<cand_pos_t> index;
    cand_pos_t total_length = ((n + 1) * (n + 2)) / 2;
    index.resize(n + 1);
    index[1] = 0;
    for (cand_pos_t i = 2; i <= n; i++) {
        index[i] = index[i - 1] + (n + 1) - i + 1;
    }
    std::vector<pf_t> M;
    std::vector<pf_t> BE;
    std::vector<pf_t> WMB;
    std::vector<pf_t> WMBP;
    std::vector<pf_t> BE_linear;
    M.resize(total_length, 0);
    BE.resize(total_length,0);
    WMB.resize(total_length,0);
    WMBP.resize(total_length,0);
    BE_linear.resize(n+1,0);

    // BlossomMatching BM(n+1,pl);
    // MEA = BM.max_weight_matching(pu,structure,tree);

    // sort by i then by j
    std::sort(pp.begin(),pp.end(),comp_plist);
    std::sort(plpk.begin(),plpk.end(),comp_plist);
    
    std::vector<cand_list_t> CL;
    CL.resize(n + 1);
    std::vector<cand_list_t> CLPK;
    CLPK.resize(n + 1);
    cand_pos_t index2 = 0;
    cand_pos_t index_PK = 0;
    cand_pos_t index_BE = 0;
    pf_t EA = 0.0;
    for(cand_pos_t i = n; i>0;--i){
        cand_pos_t ii = get_index(index,i,i);
        M[ii] = pu[i];
        if (pp[index_BE].i == i){
            BE_linear[i] = 2 * gamma * pp[index_BE].p;
            BE[ii] = 2 * gamma * pp[index_BE].p;
            index_BE++;
        }
        if ((pp[index2].i == i) && (pp[index2].i > pp[index2].j)) ++index2;
        for(cand_pos_t j = i;j<=n;++j){
            cand_pos_t ij = get_index(index,i,j);
            cand_pos_t ijm1 = get_index(index,i,j-1);
            if(tree.tree[j].pair<0) M[ij] = M[ijm1] + pu[j];

            if(tree.weakly_closed(i,j)){
                for (auto it = CL[j].begin(); CL[j].end() != it && it->first >= i; ++it) {
                    cand_pos_t k = it->first;
                    cand_pos_t ikm1 = get_index(index,i,k-1);
                    EA = it->second + M[ikm1];
                    M[ij] = std::max(M[j], EA);
                }

                for (auto it = CLPK[j].begin(); CLPK[j].end() != it && it->first >= i; ++it) {
                    cand_pos_t k = it->first;
                    cand_pos_t ikm1 = get_index(index,i,k-1);
                    EA    = it->second + M[ikm1];
                    M[ij] = std::max(M[j], EA);
                }
            }

            if ((pp[index2].i == i) && (pp[index2].j == j)) {
                cand_pos_t ipjm1 = get_index(index,i+1,j-1);
                EA = M[ipjm1];
                EA += 2 * gamma * pp[index2].p;
                if (M[ij] < EA) {
                    M[ij] = EA;
                    register_candidate(CL,i,j,EA);
                }
                ++index2;
            }
            if ((plpk[index_PK].i == i) && (plpk[index_PK].j == j)) {
                cand_pos_t bp_ij = tree.bp(i,j);
                cand_pos_t Bp_ij = tree.Bp(i,j);
                cand_pos_t b_ij = tree.b(i,j);
                cand_pos_t B_ij = tree.B(i,j);
                EA = 0;
                if(Bp_ij>=0 && B_ij>=0 && bp_ij<0){
                    cand_pos_t ip1Bpm1 = get_index(index,i+1,Bp_ij-1);
                    cand_pos_t Bp1jm1 = get_index(index,B_ij+1,j-1);
                    EA = M[ip1Bpm1] + M[Bp1jm1];
                }
                if(b_ij>=0 && bp_ij>=0 && Bp_ij<0){
                    cand_pos_t ip1bm1 = get_index(index,i+1,b_ij-1);
                    cand_pos_t bp1jm1 = get_index(index,bp_ij+1,j-1);
                    EA = M[ip1bm1] + M[bp1jm1];
                }
                if(Bp_ij>=0 && B_ij>=0 && bp_ij>0 && b_ij>0){
                    cand_pos_t ip1Bpm1 = get_index(index,i+1,Bp_ij-1);
                    cand_pos_t Bp1bm1 = get_index(index,B_ij+1,b_ij-1);
                    cand_pos_t bp1jm1 = get_index(index,bp_ij+1,j-1);
                    EA = M[ip1Bpm1] + M[Bp1bm1] + M[bp1jm1];
                }
                cand_pos_t ip1jm1 = get_index(index,i+1,j-1);
                EA = std::max(EA,M[ip1jm1]);
                EA += 2 * gamma * plpk[index_PK].p;
                if (M[ij] < EA) {
                    M[ij] = EA;
                }
                ++index_PK;
            }

            if (!((j-i-1) <= TURN || (tree.tree[i].pair >= -1 && tree.tree[i].pair > j) || (tree.tree[j].pair >= -1 && tree.tree[j].pair < i) || (tree.tree[i].pair >= -1 && tree.tree[i].pair < i ) || (tree.tree[j].pair >= -1 && j < tree.tree[j].pair))){
                // WMBP
                pf_t tmp = 0;
                // 1
                if (tree.tree[j].pair < 0){
                    cand_pos_t b_ij = tree.b(i,j);
                    for (auto it = plpk.begin(); plpk.end() != it; ++it){
                        if(j == it->j){
                            cand_pos_t l = it->i;
                            cand_pos_t bp_il = tree.bp(i,l);
			                cand_pos_t Bp_lj = tree.Bp(l,j);
                            int ext_case = compute_exterior_cases(l,j,tree);
                            if((b_ij > 0 && l < b_ij) || (b_ij<0 && ext_case == 0)){
                                if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){
                                    cand_pos_t B_lj = tree.B(l,j);
                                    if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l+TURN <=j){
                                        cand_pos_t BE_index = get_index(index,tree.tree[B_lj].pair,tree.tree[Bp_lj].pair); 
                                        cand_pos_t WMBP_index = get_index(index,i,l-1);
                                        cand_pos_t VP_index = get_index(index,l,j); // VP can still be stored in M
                                        tmp = std::max(tmp,BE[BE_index] + WMBP[WMBP_index] + M[VP_index]);
                                    }
                                }
                            }
                        }
                    }
                }
                WMBP[ij] = std::max(WMBP[ij],tmp);

                // 3
                WMBP[ij] = std::max(WMBP[ij],M[ij]);

                // 4
                tmp = 0;
                if(tree.tree[j].pair < 0 && tree.tree[i].pair >= 0){
                    for (auto it = plpk.begin(); plpk.end() != it; ++it){
                        if(j == it->j){
                            cand_pos_t l = it->i;
                            cand_pos_t bp_il = tree.bp(i,l);
                            if(bp_il >= 0 && bp_il < n && l+TURN <= j){
                                cand_pos_t index_BE = get_index(index,i,bp_il);
                                cand_pos_t index_WI = get_index(index,bp_il+1,l-1);
                                cand_pos_t index_VP = get_index(index,l,j);
                                tmp = std::max(tmp,BE[index_BE] + M[index_WI] + M[index_VP]);
                            }
                        }  
                    }
                }
                WMBP[ij] = std::max(WMBP[ij],tmp);

                // WMB
                tmp = 0;
                if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i){
                    cand_pos_t bp_j = tree.tree[j].pair;
		            for (auto it = plpk.begin(); plpk.end() != it && it->j >= bp_j; ++it){
                        cand_pos_t l = it->j;
                        cand_pos_t Bp_lj = tree.Bp(l,j);
                        cand_pos_t BE_index = get_index(index,bp_j,Bp_lj);
                        cand_pos_t WMBP_index = get_index(index,i,l);
                        cand_pos_t M_index = get_index(index,l+1,Bp_lj-1);
                        tmp = std::max(tmp,BE[BE_index] + WMBP[WMBP_index] + M[M_index]);
                    }
                }
                tmp = std::max(tmp,WMBP[ij]);
                WMB[ij] = tmp;
                
            }
            if(M[ij] < WMB[ij]){
                M[ij] = std::max(M[ij],WMB[ij]);
                register_candidate(CLPK,i,j,M[ij]);
            }
            
            cand_pos_t ip = tree.tree[i].pair; // i's pair ip should be right side so ip = )
            cand_pos_t jp = tree.tree[j].pair; // j's pair jp should be right side so jp = )
            //    ( (.     (.     ).  ).  )
            // .  i l.     j.     jp. lp. ip
            if(i>=1 && i<j && jp<ip && j<jp && i<ip && j<=n){
                //  cand_pos_t ijm1 = index[i] + j - i;
                if(i == j && ip == jp){
                    BE[ij] = BE_linear[i];
                }
                else if(i+1 == j && ip-1 == jp){
                    cand_pos_t ip1jm1 = index[i+1] + j-1 - (i+1);
                    BE[ij] = BE[ip1jm1] + BE_linear[i];
                } else{
                    pf_t m2 = 0;
                    for (cand_pos_t l = i+1; l<= j ; l++){
                        if (tree.tree[l].pair >= -1 && jp <= tree.tree[l].pair && tree.tree[l].pair < ip){
                            cand_pos_t lp = tree.tree[l].pair;
                            cand_pos_t ip1lm1 = get_index(index,i+1,l-1);
                            cand_pos_t lpp1ipm1 = get_index(index,lp+1,ip-1);
                            cand_pos_t lj = get_index(index,l,j);
                            pf_t tmp = BE_linear[i] + BE[lj];
                            if(!(i+1>l-1)) tmp += M[ip1lm1];
                            if(!(lp+1>ip-1)) tmp += M[lpp1ipm1];
                            m2 = std::max(m2,tmp);
                        }
                    }
                    BE[ij] = m2;
                }
            }
        }
    }
    cand_pos_t index_1n = get_index(index,1,n);
    MEA = M[index_1n];

    MEAdat bdat(index,pp,plpk,pu,M,BE,WMBP,gamma,CL,CLPK,structure);
    mea_backtrack(bdat,tree,1,n,0);


    this->MEA_structure = bdat.structure;

    return MEA;
}
void mea_backtrack_vp(MEAdat &bdat,sparse_tree &tree,cand_pos_t i,cand_pos_t j){
    bdat.structure[i-1] = '[';
    bdat.structure[j-1] = ']';
    cand_pos_t ij = get_index(bdat.index,i,j);
    cand_pos_t ip1jm1 = get_index(bdat.index,i+1,j-1);
    cand_pos_t bp_ij = tree.bp(i,j);
    cand_pos_t Bp_ij = tree.Bp(i,j);
    cand_pos_t b_ij = tree.b(i,j);
    cand_pos_t B_ij = tree.B(i,j);
    pf_t prec = std::numeric_limits<double>::epsilon() * bdat.M[ij];
    pf_t EA = 0.0;
    for (auto it = bdat.plpk.begin(); bdat.plpk.end() != it; ++it){
        if(i == it->i && j == it->j){
            EA = 2 * bdat.gamma * it->p;
        } 
    }
    if (bdat.M[ij] <= bdat.M[ip1jm1] + EA + prec) {
        mea_backtrack_vp(bdat,tree,i+1,j-1);
        return;
    }

    if(Bp_ij>=0 && B_ij>=0 && bp_ij<0){
        mea_backtrack(bdat,tree,i+1,Bp_ij-1,0);
        mea_backtrack(bdat,tree,B_ij+1,j-1,0);
        return;
    }
    if(b_ij>=0 && bp_ij>=0 && Bp_ij<0){
        mea_backtrack(bdat,tree,i+1,b_ij-1,0);
        mea_backtrack(bdat,tree,bp_ij+1,j-1,0);
        return;
    }
    if(Bp_ij>=0 && B_ij>=0 && bp_ij>0 && b_ij>0){
        mea_backtrack(bdat,tree,i+1,Bp_ij-1,0);
        mea_backtrack(bdat,tree,B_ij+1,b_ij-1,0);
        mea_backtrack(bdat,tree,bp_ij+1,j-1,0);
        return;
    }
}
void mea_backtrack_pk(MEAdat &bdat,sparse_tree &tree,cand_pos_t i,cand_pos_t j, pf_t e){
    // WMB
    if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i){
        cand_pos_t bp_j = tree.tree[j].pair;
        cand_pos_t Bp_lj_final = -1;
        cand_pos_t l_final = -1;
        cand_pos_t tmp = 0;
        for (auto it = bdat.plpk.begin(); bdat.plpk.end() != it && it->j >= bp_j; ++it){
            cand_pos_t l = it->j;
            cand_pos_t Bp_lj = tree.Bp(l,j);
            cand_pos_t BE_index = get_index(bdat.index,bp_j,Bp_lj);
            cand_pos_t WMBP_index = get_index(bdat.index,i,l);
            cand_pos_t M_index = get_index(bdat.index,l+1,Bp_lj-1);
            pf_t en = bdat.BE[BE_index] + bdat.WMBP[WMBP_index] + bdat.M[M_index];
            if(en>tmp){
                tmp = en;
                l_final = l;
                Bp_lj_final = Bp_lj;
            }
        }
        if(e == tmp){
            cand_pos_t WMBP_index = get_index(bdat.index,i,l_final);
            mea_backtrack(bdat,tree,l_final+1,Bp_lj_final-1,0);
            mea_backtrack_pk(bdat,tree,i,l_final,bdat.WMBP[WMBP_index]);
            for(cand_pos_t k = Bp_lj_final; k<=bp_j;++k){
                if(tree.tree[k].pair>0){
                    bdat.structure[k-1] = ')';
                    bdat.structure[tree.tree[k].pair-1] = '(';
                }
            }
        }
    }

    //WMBP 1
    if (tree.tree[j].pair < 0){
        cand_pos_t b_ij = tree.b(i,j);
        pf_t tmp = 0;
        cand_pos_t Bp_lj_final = -1;
        cand_pos_t B_lj_final = -1;
        cand_pos_t l_final = -1;
        for (auto it = bdat.plpk.begin(); bdat.plpk.end() != it; ++it){
            if(j == it->j){
                cand_pos_t l = it->i;
                cand_pos_t bp_il = tree.bp(i,l);
                cand_pos_t Bp_lj = tree.Bp(l,j);
                int ext_case = compute_exterior_cases(l,j,tree);
                if((b_ij > 0 && l < b_ij) || (b_ij<0 && ext_case == 0)){
                    if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){
                        cand_pos_t B_lj = tree.B(l,j);
                        if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l+TURN <=j){
                            cand_pos_t BE_index = get_index(bdat.index,tree.tree[B_lj].pair,tree.tree[Bp_lj].pair);
                            cand_pos_t WMBP_index = get_index(bdat.index,i,l-1);
                            cand_pos_t VP_index = get_index(bdat.index,l,j); // VP can still be stored in M
                            pf_t en = bdat.BE[BE_index] + bdat.WMBP[WMBP_index] + bdat.M[VP_index];
                            if(en>tmp){
                                tmp = en;
                                l_final = l;
                                B_lj_final = B_lj;
                                Bp_lj_final = Bp_lj;
                            }
                        }
                    }
                }
            }
        }
        if(e==tmp){
            cand_pos_t WMBP_index = get_index(bdat.index,i,l_final-1);
            mea_backtrack_pk(bdat,tree,i,l_final,bdat.WMBP[WMBP_index]);
            mea_backtrack_vp(bdat,tree,l_final,j);

            for(cand_pos_t k = Bp_lj_final; k<=B_lj_final;++k){
                if(tree.tree[k].pair>0){
                    bdat.structure[k-1] = ')';
                    bdat.structure[tree.tree[k].pair-1] = '(';
                }
            }
        }
    }
    //WMBP 2

    //WMBP 3
    cand_pos_t ij = get_index(bdat.index,i,j);
    if(!tree.weakly_closed(i,j) && e == bdat.M[ij]){
        mea_backtrack_vp(bdat,tree,i,j);
    }
    //WMBP 4
    if(tree.tree[j].pair < 0 && tree.tree[i].pair >= 0){
        pf_t tmp = 0;
        cand_pos_t bp_il_final = -1;
        cand_pos_t l_final = -1;
        for (auto it = bdat.plpk.begin(); bdat.plpk.end() != it; ++it){
            if(j == it->j){
                cand_pos_t l = it->i;
                cand_pos_t bp_il = tree.bp(i,l);
                if(bp_il >= 0 && l+TURN <= j){
                    cand_pos_t index_BE = get_index(bdat.index,i,bp_il);
                    cand_pos_t index_WI = get_index(bdat.index,bp_il+1,l-1);
                    cand_pos_t index_VP = get_index(bdat.index,l,j);
                    pf_t en = bdat.BE[index_BE] + bdat.M[index_WI] + bdat.M[index_VP];
                    if(en>tmp){
                        tmp = en;
                        bp_il_final = bp_il;
                        l_final = l;
                    }
                }
            }  
        }
        if(e==tmp){
            mea_backtrack(bdat,tree,bp_il_final+1,l_final-1,0);
            mea_backtrack_vp(bdat,tree,l_final,j);
            for(cand_pos_t k = i; k<=bp_il_final;++k){
                if(tree.tree[k].pair>0){
                    bdat.structure[k-1] = '(';
                    bdat.structure[tree.tree[k].pair-1] = ')';
                }
            }
            return;
        }
    }
}

void mea_backtrack(MEAdat &bdat,sparse_tree &tree,cand_pos_t i,cand_pos_t j, int pair){

    int fail  = 1;

    if(tree.tree[i].pair == j){
        bdat.structure[i - 1]  = '(';
        bdat.structure[j - 1]  = ')';
        ++i;
        --j;
    }
    if (pair) {
        /*
        * if pair == 1, insert pair 
        */
        bdat.structure[i - 1]  = '(';
        bdat.structure[j - 1]  = ')';
        ++i;
        --j;
    }
    cand_pos_t ij = get_index(bdat.index,i,j);
    cand_pos_t ijm1 = get_index(bdat.index,i,j-1);

    pf_t prec = std::numeric_limits<double>::epsilon() * bdat.M[ij];
    /* Mi values are filled, do the backtrace */
    while ((j > i) && (bdat.M[ij] <= (bdat.M[ijm1] + bdat.pu[j] + prec))) {
        bdat.structure[j - 1] = '.';
        --j;
        --ij; // Since ij is index[i] + j - i, if I subtract 1, it makes it j-1
        --ijm1;
    }

    for (auto it = bdat.CL[j].begin(); bdat.CL[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        cand_pos_t ikm1 = get_index(bdat.index,i,k-1);
        if (bdat.M[ij] <= it->second + bdat.M[ikm1] + prec) {
            if (k > i + 3) mea_backtrack(bdat,tree, i, k-1,0);

            mea_backtrack(bdat,tree, k, j,1);
            fail = 0;
        }
    }
    for (auto it = bdat.CLPK[j].begin(); bdat.CLPK[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        cand_pos_t ikm1 = get_index(bdat.index,i,k-1);
        if (bdat.M[ij] <= it->second + bdat.M[ikm1] + prec) {
            if (k > i + 3) mea_backtrack(bdat,tree, i, k-1,0);

            mea_backtrack_pk(bdat,tree, k, j,it->second);
            fail = 0;
        }
    }
    if (fail && j > i){
        printf("backtrack failed for MEA at %d and %d\n",i,j);
        exit(0);
    }
}