#ifndef PART_FUNC
#define PART_FUNC
#include "base_types.hh"
#include "sparse_tree.hh"
#include <cstring>
#include <string>
#include <vector>

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
#include "ViennaRNA/params/default.h"
}


class W_final_pf{

    public:
        W_final_pf(std::string seq,bool pk_only, int dangle, double energy);
        // constructor for the restricted mfe case

        ~W_final_pf ();
        // The destructor

        double hfold_pf (sparse_tree &tree);

        vrna_exp_param_t *exp_params_;

        pf_t get_energy (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return V[ij]; }
        pf_t get_energy_WM (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WM[ij]; }
        pf_t get_energy_WMv (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WMv[ij]; }
        pf_t get_energy_WMp (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WMp[ij]; }

        pf_t get_energy_WI (cand_pos_t i, cand_pos_t j) { if (i>j) return 1; cand_pos_t ij = index[i]+j-i; return WI[ij]; }
        pf_t get_energy_WIP (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WIP[ij]; }
        pf_t get_energy_VP (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return VP[ij]; }
        pf_t get_energy_VPL (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return VPL[ij]; }
        pf_t get_energy_VPR (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return VPR[ij]; }
        pf_t get_energy_WMB (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WMB[ij]; }
        pf_t get_energy_WMBP (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WMBP[ij]; }
        pf_t get_energy_WMBW (cand_pos_t i, cand_pos_t j) { if (i>=j) return 0; cand_pos_t ij = index[i]+j-i; return WMBW[ij]; }
        pf_t get_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree){
        // Hosna, March 16, 2012,
        // i and j should be at least 3 bases apart
            if (j-i>= TURN && i >= 1 && i <= ip && ip < jp && jp <= j && j <=n && tree.tree[i].pair >=0 && tree.tree[j].pair >= 0 && tree.tree[ip].pair >= 0 && tree.tree[jp].pair >= 0 && tree.tree[i].pair == j && tree.tree[j].pair == i && tree.tree[ip].pair == jp && tree.tree[jp].pair == ip){
                // if(i == ip && j == jp && i<j){
                //     return 1;
                // }
                cand_pos_t iip = index[i]+ip-i;

                return BE[iip];
            }else{
                return 0;
            }
        }
        

    private:
        std::string seq;
        bool pk_free;
        cand_pos_t n;
        std::vector<cand_pos_t> index;

        short *S_;
        short *S1_;

        std::vector<pf_t> V;
        std::vector<pf_t> WMv;
        std::vector<pf_t> WMp;
        std::vector<pf_t> WM;
        std::vector<pf_t> W;

        std::vector<pf_t> WI;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot)
        std::vector<pf_t> VP;				// the loop corresponding to the pseudoknotted region of WMB
        std::vector<pf_t> VPL;				// the loop corresponding to the pseudoknotted region of WMB
        std::vector<pf_t> VPR;				// the loop corresponding to the pseudoknotted region of WMB
        std::vector<pf_t> WMB;				// the main loop for pseudoloops and bands
        std::vector<pf_t> WMBP; 				// the main loop to calculate WMB
        std::vector<pf_t> WMBW;
        std::vector<pf_t> WIP;				// the loop corresponding to WI'
        std::vector<pf_t> BE;				// the loop corresponding to BE

        std::vector<pf_t> scale;
        std::vector<pf_t> expMLbase;
        std::vector<pf_t> expcp_pen;
        std::vector<pf_t> expPUP_pen;

        void rescale_pk_globals();

        void exp_params_rescale(double mfe);

        void compute_energy_restricted(cand_pos_t i,cand_pos_t j,sparse_tree &tree);

        void compute_WMv_WMp(cand_pos_t i, cand_pos_t j, std::vector<Node> &tree);

        void compute_energy_WM_restricted (cand_pos_t i, cand_pos_t j, sparse_tree &tree);

        void compute_pk_energies(cand_pos_t i,cand_pos_t j,sparse_tree &tree);

        void compute_WI(cand_pos_t i,cand_pos_t j,sparse_tree &tree);

        void compute_WIP(cand_pos_t i,cand_pos_t j,sparse_tree &tree);

        void compute_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

        void compute_VPL(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

        void compute_VPR(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

        void compute_WMBW(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

        void compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

        void compute_WMB(cand_pos_t  i, cand_pos_t  j, sparse_tree &tree);

        void compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree);


        pf_t exp_Extloop(cand_pos_t i, cand_pos_t j);

        pf_t exp_MLstem(cand_pos_t i, cand_pos_t j);

        pf_t exp_Mbloop(cand_pos_t i, cand_pos_t j);

        pf_t HairpinE(cand_pos_t i, cand_pos_t j);

        pf_t compute_internal_restricted(cand_pos_t i, cand_pos_t j,std::vector<int> &up);

        pf_t compute_energy_VM_restricted (cand_pos_t i, cand_pos_t j, std::vector<Node> &tree);

        pf_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);

        pf_t get_e_stP(cand_pos_t i, cand_pos_t j);

        pf_t get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j);

        int compute_exterior_cases(cand_pos_t l, cand_pos_t j, const sparse_tree &tree);
};

#endif