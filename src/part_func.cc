#include "part_func.hh"
#include "dot_plot.hh"
#include "h_externs.hh"
#include "pf_globals.hh"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>

#define debug 0
/*
 * If the global use_mfelike_energies flag is set, truncate doubles to int
 * values and cast back to double. This makes the energy parameters of the
 * partition (folding get_scaled_exp_params()) compatible with the mfe folding
 * parameters (get_scaled_exp_params()), e.g. for explicit partition function
 * computations.
 */
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))
/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT) ((dH) - ((dH) - (dG)) * dT)

/*
 * Rescale Free energy contribution according to deviation of temperature from measurement conditions
 * and convert it to Boltzmann Factor for specific kT
 */
#define RESCALE_BF(dG, dH, dT, kT) (exp(-TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) * 10. / kT))

W_final_pf::W_final_pf(std::string &seq, std::string &MFE_structure, bool pk_free, int dangle, double energy, int num_samples, bool PSplot)
    : exp_params_(scale_pf_parameters()) {
    this->seq = seq;
    this->MFE_structure = MFE_structure;
    this->n = seq.length();
    this->pk_free = pk_free;
    this->PSplot = PSplot;
    this->num_samples = num_samples;

    make_pair_matrix();
    exp_params_->model_details.dangles = dangle;
    S_ = encode_sequence(seq.c_str(), 0);
    S1_ = encode_sequence(seq.c_str(), 1);

    index.resize(n + 1);
    scale.resize(n + 1);
    expMLbase.resize(n + 1);
    expcp_pen.resize(n + 1);
    expPUP_pen.resize(n + 1);
    cand_pos_t total_length = ((n + 1) * (n + 2)) / 2;
    index[1] = 0;
    for (cand_pos_t i = 2; i <= n; i++)
        index[i] = index[i - 1] + (n + 1) - i + 1;
    // Allocate space
    V.resize(total_length, 0);
    VM.resize(total_length, 0);
    WM.resize(total_length, 0);
    WMv.resize(total_length, 0);
    WMp.resize(total_length, 0);
    // W.resize(n+1,1);

    // PK
    WIP.resize(total_length, 0);
    VP.resize(total_length, 0);
    VPL.resize(total_length, 0);
    VPR.resize(total_length, 0);
    WMB.resize(total_length, 0);
    WMBP.resize(total_length, 0);
    WMBW.resize(total_length, 0);
    BE.resize(total_length, 0);

    rescale_pk_globals();
    exp_params_rescale(energy);
    W.resize(n + 1, scale[1]);
    WI.resize(total_length, scale[1]);
}

W_final_pf::~W_final_pf() {}

void W_final_pf::exp_params_rescale(double mfe) {
    double e_per_nt, kT;
    kT = exp_params_->kT;

    e_per_nt = mfe * 1000. / this->n;

    exp_params_->pf_scale = exp(-(exp_params_->model_details.sfact * e_per_nt) / kT);

    if (exp_params_->pf_scale < 1.) exp_params_->pf_scale = 1.;

    exp_params_->pf_scale = 1.;

    this->scale[0] = 1.;
    this->scale[1] = (pf_t)(1. / exp_params_->pf_scale);
    this->expMLbase[0] = 1;
    this->expMLbase[1] = (pf_t)(exp_params_->expMLbase / exp_params_->pf_scale);

    this->expcp_pen[0] = 1;
    this->expcp_pen[1] = (pf_t)(expcp_penalty / exp_params_->pf_scale);
    this->expPUP_pen[0] = 1;
    this->expPUP_pen[1] = (pf_t)(expPUP_penalty / exp_params_->pf_scale);

    for (cand_pos_t i = 2; i <= this->n; i++) {
        this->scale[i] = this->scale[i / 2] * this->scale[i - (i / 2)];
        this->expMLbase[i] = (pf_t)pow(exp_params_->expMLbase, (double)i) * this->scale[i];
        this->expcp_pen[i] = (pf_t)pow(expcp_penalty, (double)i) * this->scale[i];
        this->expPUP_pen[i] = (pf_t)pow(expPUP_penalty, (double)i) * this->scale[i];
    }
}

void W_final_pf::rescale_pk_globals() {
    double kT = exp_params_->model_details.betaScale * (exp_params_->model_details.temperature + K0) * GASCONST; /* kT in cal/mol  */
    double TT = (exp_params_->model_details.temperature + K0) / (Tmeasure);
    int pf_smooth = exp_params_->model_details.pf_smooth;

    expPS_penalty = RESCALE_BF(PS_penalty, PS_penalty * 3, TT, kT);
    expPSM_penalty = RESCALE_BF(PSM_penalty, PSM_penalty * 3, TT, kT);
    expPSP_penalty = RESCALE_BF(PSP_penalty, PSP_penalty * 3, TT, kT);
    expPB_penalty = RESCALE_BF(PB_penalty, PB_penalty * 3, TT, kT);
    expPUP_penalty = RESCALE_BF(PUP_penalty, PUP_penalty * 3, TT, kT);
    expPPS_penalty = RESCALE_BF(PPS_penalty, PPS_penalty * 3, TT, kT);

    expa_penalty = RESCALE_BF(a_penalty, ML_closingdH, TT, kT);
    expb_penalty = RESCALE_BF(b_penalty, ML_interndH, TT, kT);
    expc_penalty = RESCALE_BF(c_penalty, ML_BASEdH, TT, kT);

    expap_penalty = RESCALE_BF(ap_penalty, ap_penalty * 3, TT, kT);
    expbp_penalty = RESCALE_BF(bp_penalty, bp_penalty * 3, TT, kT);
    expcp_penalty = RESCALE_BF(cp_penalty, cp_penalty * 3, TT, kT);
}

/**
 * In cases where the band border is not found, if specific cases are met, the value is Inf(i.e n) not -1.
 * When applied to WMBP, if all cases are 0, then we can proceed with WMBP
 * Mateo Jan 2025: Added to Fix WMBP problem
 */
int W_final_pf::compute_exterior_cases(cand_pos_t l, cand_pos_t j, sparse_tree &tree) {
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

inline double W_final_pf::to_Energy(pf_t energy, cand_pos_t length) {
    return ((-log(energy) - length * log(exp_params_->pf_scale)) * exp_params_->kT / 1000.0);
}

double W_final_pf::hfold_pf(sparse_tree &tree) {

    for (cand_pos_t i = n; i >= 1; --i) {
        for (cand_pos_t j = i; j <= n; ++j) {

            const bool evaluate = tree.weakly_closed(i, j);
            const pair_type ptype_closing = pair[S_[i]][S_[j]];
            const bool restricted = tree.tree[i].pair == -1 || tree.tree[j].pair == -1;

            if (ptype_closing > 0 && evaluate && !restricted) compute_energy_restricted(i, j, tree);

            if (!pk_free) compute_pk_energies(i, j, tree);

            compute_WMv_WMp(i, j, tree.tree);
            compute_energy_WM_restricted(i, j, tree);
        }
    }
    for (cand_pos_t j = TURN + 1; j <= n; j++) {
        pf_t contributions = 0;
        if (tree.tree[j].pair < 0) contributions += W[j - 1] * scale[1];
        if (tree.weakly_closed(1, j)) {
            for (cand_pos_t k = 1; k <= j - TURN - 1; ++k) {
                if (tree.weakly_closed(1, k - 1)) {
                    pf_t acc = (k > 1) ? W[k - 1] : 1; // keep as 0 or 1?
                    contributions += acc * get_energy(k, j) * exp_Extloop(k, j);
                    if (k == 1 || tree.weakly_closed(k, j)) contributions += acc * get_energy_WMB(k, j) * expPS_penalty;
                }
            }
        }
        W[j] = contributions;
    }

    pf_t energy = to_Energy(W[n], n);

    // Base pair probability
    structure = std::string(n, '.');
    std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples;
    std::unordered_map<std::string, int> structures;
    for (cand_pos_t i = 0; i < num_samples; ++i) {
        std::string structure(n, '.');
        Sample_W(1, n, structure, samples, tree);
        structures[structure]++;
    }
    // for (auto const& [key, val] : structures){
	// 	std::cout << key << ':' << val << std::endl;
    // }
    pairing_tendency(samples, tree);
    frequency = (pf_t)structures[MFE_structure] / num_samples;
    std::cout << "frequency of MFE structure in ensemble: " << frequency << std::endl;

    if (PSplot) {
        create_dot_plot(seq, tree.tree, MFE_structure, samples, num_samples);
    }

    return energy;
}

pf_t W_final_pf::exp_Extloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];

    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i - 1] : -1;
        base_type sj1 = j < n ? S_[j + 1] : -1;
        return exp_E_ExtLoop(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_ExtLoop(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::exp_MLstem(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];
    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i - 1] : -1;
        base_type sj1 = j < n ? S_[j + 1] : -1;
        return exp_E_MLstem(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::exp_Mbloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[j]][S_[i]];
    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i + 1] : -1;
        base_type sj1 = j < n ? S_[j - 1] : -1;
        return exp_E_MLstem(tt, sj1, si1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::HairpinE(cand_pos_t i, cand_pos_t j) {

    const int ptype_closing = pair[S_[i]][S_[j]];
    if (ptype_closing == 0) return 0;
    pf_t e_h = static_cast<pf_t>(exp_E_Hairpin(j - i - 1, ptype_closing, S1_[i + 1], S1_[j - 1], &seq.c_str()[i - 1], exp_params_));
    e_h *= scale[j - i + 1];
    return e_h;
}

pf_t W_final_pf::compute_internal_restricted(cand_pos_t i, cand_pos_t j, std::vector<int> &up) {
    pf_t v_iloop = 0;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1);
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (cand_pos_t k = i + 1; k <= max_k; ++k) {
        if ((up[k - 1] >= (k - i - 1))) {
            cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
            for (cand_pos_t l = j - 1; l >= min_l; --l) {
                if (up[j - 1] >= (j - l - 1)) {
                    pf_t v_iloop_kl = get_energy(k, l)
                                      * exp_E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1],
                                                      S1_[k - 1], S1_[l + 1], exp_params_);
                    cand_pos_t u1 = k - i - 1;
                    cand_pos_t u2 = j - l - 1;
                    v_iloop_kl *= scale[u1 + u2 + 2];
                    v_iloop += v_iloop_kl;
                }
            }
        }
    }

    return v_iloop;
}

void W_final_pf::compute_WMv_WMp(cand_pos_t i, cand_pos_t j, std::vector<Node> &tree) {
    if (j - i - 1 < TURN) return;
    cand_pos_t ij = index[(i)] + (j) - (i);

    pf_t WMv_contributions = 0;
    pf_t WMp_contributions = 0;

    WMv_contributions += (get_energy(i, j) * exp_MLstem(i, j));
    WMp_contributions += (get_energy_WMB(i, j) * expPSM_penalty * expb_penalty);
    if (tree[j].pair < 0) {
        WMv_contributions += (get_energy_WMv(i, j - 1) * expMLbase[1]);
        WMp_contributions += (get_energy_WMp(i, j - 1) * expMLbase[1]);
    }
    WMv[ij] = WMv_contributions;
    WMp[ij] = WMp_contributions;
}

void W_final_pf::compute_energy_WM_restricted(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    if (j - i + 1 < 4) return;
    pf_t contributions = 0;
    cand_pos_t ij = index[(i)] + (j) - (i);
    cand_pos_t ijminus1 = index[(i)] + (j)-1 - (i);

    for (cand_pos_t k = i; k < j - TURN; ++k) {
        pf_t qbt1 = get_energy(k, j) * exp_MLstem(k, j);
        pf_t qbt2 = get_energy_WMB(k, j) * expPSM_penalty * expb_penalty;
        bool can_pair = tree.up[k - 1] >= (k - i);
        if (can_pair) contributions += (static_cast<pf_t>(expMLbase[k - i]) * qbt1);
        if (can_pair) contributions += (static_cast<pf_t>(expMLbase[k - i]) * qbt2);
        contributions += (get_energy_WM(i, k - 1) * qbt1);
        contributions += (get_energy_WM(i, k - 1) * qbt2);
    }
    if (tree.tree[j].pair < 0) contributions += WM[ijminus1] * expMLbase[1];
    WM[ij] = contributions;
}

pf_t W_final_pf::compute_energy_VM_restricted(cand_pos_t i, cand_pos_t j, std::vector<int> &up) {
    pf_t contributions = 0;
    cand_pos_t ij = index[(i)] + (j) - (i);
    for (cand_pos_t k = i + 1; k <= j - TURN - 1; ++k) {
        contributions += (get_energy_WM(i + 1, k - 1) * get_energy_WMv(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        contributions += (get_energy_WM(i + 1, k - 1) * get_energy_WMp(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        if (up[k - 1] >= (k - (i + 1)))
            contributions += (expMLbase[k - i - 1] * get_energy_WMp(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
    }

    contributions *= scale[2];
    VM[ij] = contributions;
    return contributions;
}

void W_final_pf::compute_energy_restricted(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {

    cand_pos_t ij = index[i] + j - i;

    const bool unpaired = (tree.tree[i].pair < -1 && tree.tree[j].pair < -1);
    const bool paired = (tree.tree[i].pair == j && tree.tree[j].pair == i);

    pf_t contributions = 0;

    if (paired || unpaired) // if i and j can pair
    {
        bool canH = !(tree.up[j - 1] < (j - i - 1));
        if (canH) contributions += HairpinE(i, j);

        contributions += compute_internal_restricted(i, j, tree.up);

        contributions += compute_energy_VM_restricted(i, j, tree.up);
    }
    V[ij] = contributions;
}

void W_final_pf::compute_pk_energies(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {

    cand_pos_t ij = index[i] + j - i;
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    bool weakly_closed_ij = tree.weakly_closed(i, j);

    if ((i == j || j - i < 4 || weakly_closed_ij)) {
        VP[ij] = 0;
        VPL[ij] = 0;
        VPR[ij] = 0;
    } else {
        if (ptype_closing > 0 && tree.tree[i].pair < -1 && tree.tree[j].pair < -1) compute_VP(i, j, tree);
        if (tree.tree[j].pair < -1) compute_VPL(i, j, tree);
        if (tree.tree[j].pair < j) compute_VPR(i, j, tree);
    }

    if (!((j - i - 1) <= TURN || (tree.tree[i].pair >= -1 && tree.tree[i].pair > j) || (tree.tree[j].pair >= -1 && tree.tree[j].pair < i)
          || (tree.tree[i].pair >= -1 && tree.tree[i].pair < i) || (tree.tree[j].pair >= -1 && j < tree.tree[j].pair))) {
        compute_WMBW(i, j, tree);
        compute_WMBP(i, j, tree);
        compute_WMB(i, j, tree);
    }

    if (!weakly_closed_ij) {
        WI[ij] = 0;
        WIP[ij] = 0;
    } else {
        compute_WI(i, j, tree);
        compute_WIP(i, j, tree);
    }
    cand_pos_t ip = tree.tree[i].pair; // i's pair ip should be right side so ip = )
    cand_pos_t jp = tree.tree[j].pair; // j's pair jp should be left side so jp = (
    compute_BE(i, ip, jp, j, tree);
}

void W_final_pf::compute_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {

    cand_pos_t ij = index[i] + j - i;
    pf_t contributions = 0;
    if (i == j) {
        WI[ij] = expPUP_pen[1];
        return;
    }
    // contributions += (get_energy(i,j)*expPPS_penalty);
    // contributions += (get_energy_WMB(i,j)*expPSP_penalty*expPPS_penalty);
    // If k-1<i, WI(i,k-1) = 1; if k-1==i, WI(i,k-1) = PUP_pen
    for (cand_pos_t k = i; k <= j - TURN - 1; ++k) {
        contributions += (get_energy_WI(i, k - 1) * get_energy(k, j) * expPPS_penalty);
        contributions += (get_energy_WI(i, k - 1) * get_energy_WMB(k, j) * expPSP_penalty * expPPS_penalty);
    }
    if (tree.tree[j].pair < 0) contributions += (get_energy_WI(i, j - 1) * expPUP_pen[1]);

    WI[ij] = contributions;
}

void W_final_pf::compute_WIP(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {

    cand_pos_t ij = index[i] + j - i;
    pf_t contributions = 0;
    contributions += get_energy(i, j) * expbp_penalty;
    contributions += get_energy_WMB(i, j) * expbp_penalty * expPSM_penalty;
    for (cand_pos_t k = i + 1; k < j - TURN - 1; ++k) {
        bool can_pair = tree.up[k - 1] >= (k - i);

        contributions += (get_energy_WIP(i, k - 1) * get_energy(k, j) * expbp_penalty);
        contributions += (get_energy_WIP(i, k - 1) * get_energy_WMB(k, j) * expbp_penalty * expPSM_penalty);
        if (can_pair) contributions += (expcp_pen[k - i] * get_energy(k, j) * expbp_penalty);
        if (can_pair) contributions += (expcp_pen[k - i] * get_energy_WMB(k, j) * expbp_penalty * expPSM_penalty);
    }
    if (tree.tree[j].pair < 0) contributions += (get_energy_WIP(i, j - 1) * expcp_pen[1]);
    WIP[ij] = contributions;
}

void W_final_pf::compute_VPL(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {

    cand_pos_t ij = index[i] + j - i;
    pf_t contributions = 0;

    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree.b(i, j), (cand_pos_tu)tree.Bp(i, j));
    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        bool can_pair = tree.up[k - 1] >= (k - i);
        if (can_pair) contributions += (expcp_pen[k - i] * get_energy_VP(k, j));
    }
    VPL[ij] = contributions;
}

void W_final_pf::compute_VPR(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {

    cand_pos_t ij = index[i] + j - i;
    pf_t contributions = 0;
    cand_pos_t max_i_bp = std::max(tree.B(i, j), tree.bp(i, j));
    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        bool can_pair = tree.up[j - 1] >= (j - k);
        contributions += (get_energy_VP(i, k) * get_energy_WIP(k + 1, j));
        if (can_pair) contributions += (get_energy_VP(i, k) * expcp_pen[k - i]);
    }
    VPR[ij] = contributions;
}

void W_final_pf::compute_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    cand_pos_t ij = index[i] + j - i;

    pf_t contributions = 0;

    // Borders -- added one to i and j to make it fit current bounds but also subtracted 1 from answer as the tree bounds are shifted as well
    cand_pos_t Bp_ij = tree.Bp(i, j);
    cand_pos_t B_ij = tree.B(i, j);
    cand_pos_t b_ij = tree.b(i, j);
    cand_pos_t bp_ij = tree.bp(i, j);

    if ((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) < (tree.tree[i].parent->index) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        pf_t m1 = (get_energy_WI(i + 1, Bp_ij - 1) * get_energy_WI(B_ij + 1, j - 1));
        m1 *= scale[2];
        contributions += m1;
    }

    if ((tree.tree[i].parent->index) < (tree.tree[j].parent->index) && (tree.tree[j].parent->index) > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        pf_t m2 = (get_energy_WI(i + 1, b_ij - 1) * get_energy_WI(bp_ij + 1, j - 1));
        m2 *= scale[2];
        contributions += m2;
    }

    if ((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        pf_t m3 = (get_energy_WI(i + 1, Bp_ij - 1) * get_energy_WI(B_ij + 1, b_ij - 1) * get_energy_WI(bp_ij + 1, j - 1));
        m3 *= scale[2];
        contributions += m3;
    }

    pair_type ptype_closingip1jm1 = pair[S_[i + 1]][S_[j - 1]];
    if ((tree.tree[i + 1].pair) < -1 && (tree.tree[j - 1].pair) < -1 && ptype_closingip1jm1 > 0) {
        pf_t vp_stp = (get_e_stP(i, j) * get_energy_VP(i + 1, j - 1));
        vp_stp *= scale[2];
        contributions += vp_stp;
    }

    cand_pos_t min_borders = std::min((cand_pos_tu)Bp_ij, (cand_pos_tu)b_ij);
    cand_pos_t edge_i = std::min(i + MAXLOOP + 1, j - TURN - 1);
    min_borders = std::min(min_borders, edge_i);
    for (cand_pos_t k = i + 1; k < min_borders; ++k) {
        if (tree.tree[k].pair < -1 && (tree.up[(k)-1] >= ((k) - (i)-1))) {
            cand_pos_t max_borders = std::max(bp_ij, B_ij) + 1;
            cand_pos_t edge_j = k + j - i - MAXLOOP - 2;
            max_borders = std::max(max_borders, edge_j);
            for (cand_pos_t l = j - 1; l > max_borders; --l) {
                pair_type ptype_closingkj = pair[S_[k]][S_[l]];
                if (k == i + 1 && l == j - 1) continue; // I have to add or else it will add a stP version and an eintP version to the sum
                if (tree.tree[l].pair < -1 && ptype_closingkj > 0 && (tree.up[(j)-1] >= ((j) - (l)-1))) {
                    pf_t vp_iloop_kl = (get_e_intP(i, k, l, j) * get_energy_VP(k, l));
                    cand_pos_t u1 = k - i - 1;
                    cand_pos_t u2 = j - l - 1;
                    vp_iloop_kl *= scale[u1 + u2 + 2];
                    contributions += vp_iloop_kl;
                }
            }
        }
    }

    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree.b(i, j), (cand_pos_tu)tree.Bp(i, j));
    cand_pos_t max_i_bp = std::max(tree.B(i, j), tree.bp(i, j));

    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        pf_t m6 = (get_energy_WIP(i + 1, k - 1) * get_energy_VP(k, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        m6 *= scale[2];
        contributions += m6;
    }

    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        pf_t m7 = (get_energy_VP(i + 1, k) * get_energy_WIP(k + 1, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        m7 *= scale[2];
        contributions += m7;
    }

    for (cand_pos_t k = i + 1; k < min_Bp_j; ++k) {
        pf_t m8 = (get_energy_WIP(i + 1, k - 1) * get_energy_VPR(k, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        m8 *= scale[2];
        contributions += m8;
    }

    for (cand_pos_t k = max_i_bp + 1; k < j; ++k) {
        pf_t m9 = (get_energy_VPL(i + 1, k) * get_energy_WIP(k + 1, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        m9 *= scale[2];
        contributions += m9;
    }

    VP[ij] = contributions;
}

pf_t W_final_pf::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) {
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return exp_E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], exp_params_);
}

pf_t W_final_pf::get_e_stP(cand_pos_t i, cand_pos_t j) {
    if (i + 1 == j - 1) { // TODO: do I need something like that or stack is taking care of this?
        return 0;
    }
    pf_t e_st = compute_int(i, j, i + 1, j - 1);

    return pow(e_st, e_stP_penalty);
}

pf_t W_final_pf::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j) {
    if (ip == i + 1 && jp == j - 1) return 0;
    pf_t e_int = compute_int(i, j, ip, jp);

    return pow(e_int, e_intP_penalty);
}

void W_final_pf::compute_WMBW(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    cand_pos_t ij = index[i] + j - i;

    pf_t contributions = 0;

    if (tree.tree[j].pair < j) {
        for (cand_pos_t l = i + 1; l < j; l++) {
            if (tree.tree[l].pair < 0 && tree.tree[l].parent->index > -1 && tree.tree[j].parent->index > -1
                && tree.tree[j].parent->index == tree.tree[l].parent->index) {
                contributions += get_energy_WMBP(i, l) * get_energy_WI(l + 1, j);
            }
        }
    }
    WMBW[ij] = contributions;
}

void W_final_pf::compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    cand_pos_t ij = index[i] + j - i;
    pf_t contributions = 0;

    if (tree.tree[j].pair < 0) {
        cand_pos_t b_ij = tree.b(i, j);
        for (cand_pos_t l = i + 1; l < j; ++l) {
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j, tree);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                cand_pos_t bp_il = tree.bp(i, l);
                cand_pos_t Bp_lj = tree.Bp(l, j);
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) {
                    cand_pos_t B_lj = tree.B(l, j);
                    if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l + TURN <= j) {
                        pf_t m1 = get_BE(tree.tree[B_lj].pair, B_lj, tree.tree[Bp_lj].pair, Bp_lj, tree) * get_energy_WMBP(i, l - 1)
                                  * get_energy_VP(l, j) * pow(expPB_penalty, 2);
                        contributions += m1;
                    }
                }
            }
        }
    }

    if (tree.tree[j].pair < 0) {
        cand_pos_t b_ij = tree.b(i, j);
        for (cand_pos_t l = i + 1; l < j; l++) {
            cand_pos_t bp_il = tree.bp(i, l);
            cand_pos_t Bp_lj = tree.Bp(l, j);
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j, tree);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) {
                    cand_pos_t B_lj = tree.B(l, j);
                    if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l + TURN <= j) {
                        pf_t m2 = get_BE(tree.tree[B_lj].pair, B_lj, tree.tree[Bp_lj].pair, Bp_lj, tree) * get_energy_WMBW(i, l - 1)
                                  * get_energy_VP(l, j) * pow(expPB_penalty, 2);
                        contributions += m2;
                    }
                }
            }
        }
    }

    pf_t m3 = get_energy_VP(i, j) * expPB_penalty;
    contributions += m3; // Make sure not to use non-Partition values

    if (tree.tree[j].pair < 0 && tree.tree[i].pair >= 0) {
        for (cand_pos_t l = i + 1; l < j; l++) {
            cand_pos_t bp_il = tree.bp(i, l);
            if (bp_il >= 0 && bp_il < n && l + TURN <= j) {
                if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l + TURN <= j) {
                    pf_t m4 = get_BE(i, tree.tree[i].pair, bp_il, tree.tree[bp_il].pair, tree) * get_energy_WI(bp_il + 1, l - 1) * get_energy_VP(l, j)
                              * pow(expPB_penalty, 2);
                    contributions += m4;
                }
            }
        }
    }

    WMBP[ij] = contributions;
}

void W_final_pf::compute_WMB(cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    cand_pos_t ij = index[i] + j - i;
    pf_t contributions = 0;
    // base case
    if (i == j) {
        WMB[ij] = 0;
        return;
    }

    if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i) {
        cand_pos_t bp_j = tree.tree[j].pair;
        for (cand_pos_t l = (bp_j + 1); (l < j); ++l) {
            // if(tree.tree[l].pair>0) continue;
            cand_pos_t Bp_lj = tree.Bp(l, j);
            if (Bp_lj >= 0 && Bp_lj < n) {
                contributions +=
                    get_BE(bp_j, j, tree.tree[Bp_lj].pair, Bp_lj, tree) * get_energy_WMBP(i, l) * get_energy_WI(l + 1, Bp_lj - 1) * expPB_penalty;
            }
        }
    }

    contributions += get_energy_WMBP(i, j);

    WMB[ij] = contributions;
}

void W_final_pf::compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree) {

    if (!(i >= 1 && i <= ip && ip < jp && jp <= j && j <= n && tree.tree[i].pair > 0 && tree.tree[j].pair > 0 && tree.tree[ip].pair > 0
          && tree.tree[jp].pair > 0 && tree.tree[i].pair == j && tree.tree[j].pair == i && tree.tree[ip].pair == jp
          && tree.tree[jp].pair == ip)) { // impossible cases
        return;
    }
    // (   (    (   )    )   ) //
    // i   l    ip  jp   lp  j //
    cand_pos_t iip = index[i] + ip - i;
    pf_t contributions = 0;
    // base case: i.j and ip.jp must be in G
    if (tree.tree[i].pair != j || tree.tree[ip].pair != jp) {
        BE[iip] = 0;
        return;
    }

    // base case:
    if (i == ip && j == jp && i < j) {

        BE[iip] = scale[2];
        return;
    }

    if (tree.tree[i + 1].pair == j - 1) {
        pf_t be_estp = get_e_stP(i, j) * get_BE(i + 1, j - 1, ip, jp, tree);
        be_estp *= scale[2];
        contributions += be_estp;
    }

    for (cand_pos_t l = i + 1; l <= ip; l++) {
        if (tree.tree[l].pair >= -1 && jp <= tree.tree[l].pair && tree.tree[l].pair < j) {

            cand_pos_t lp = tree.tree[l].pair;

            bool empty_region_il = (tree.up[(l)-1] >= l - i - 1);       // empty between i+1 and l-1
            bool empty_region_lpj = (tree.up[(j)-1] >= j - lp - 1);     // empty between lp+1 and j-1
            bool weakly_closed_il = tree.weakly_closed(i + 1, l - 1);   // weakly closed between i+1 and l-1
            bool weakly_closed_lpj = tree.weakly_closed(lp + 1, j - 1); // weakly closed between lp+1 and j-1

            if (empty_region_il && empty_region_lpj) { //&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
                pf_t eintp = get_e_intP(i, l, lp, j) * get_BE(l, lp, ip, jp, tree);
                cand_pos_t u1 = l - i - 1;
                cand_pos_t u2 = j - lp - 1;
                eintp *= scale[u1 + u2 + 2];
                contributions += eintp; // Added to e_intP that l != i+1 and lp != j-1 at the same time
            }
            if (weakly_closed_il && weakly_closed_lpj) {
                pf_t m3 = get_energy_WIP(i + 1, l - 1) * get_BE(l, lp, ip, jp, tree) * get_energy_WIP(lp + 1, j - 1) * expap_penalty
                          * pow(expbp_penalty, 2);
                m3 *= scale[2];
                contributions += m3;
            }
            if (weakly_closed_il && empty_region_lpj) {
                pf_t m4 = get_energy_WIP(i + 1, l - 1) * get_BE(l, lp, ip, jp, tree) * expcp_pen[j - lp - 1] * expap_penalty * pow(bp_penalty, 2);
                m4 *= scale[2];
                contributions += m4;
            }
            if (empty_region_il && weakly_closed_lpj) {
                pf_t m5 = expcp_pen[l - i - 1] * get_BE(l, lp, ip, jp, tree) * get_energy_WIP(lp + 1, j - 1) * expap_penalty * pow(bp_penalty, 2);
                m5 *= scale[2];
                contributions += m5;
            }
        }
    }

    BE[iip] = contributions;
}

/*                                BPP                                            */

inline cand_pos_t boustrophedon_at(cand_pos_t start, cand_pos_t end, cand_pos_t pos) {
    cand_pos_t count = pos - 1;
    cand_pos_t advance = (cand_pos_t)(count / 2);

    return start + (end - start) * (count % 2) + advance - (2 * (count % 2)) * advance;
}

std::vector<cand_pos_t> boustrophedon(cand_pos_t start, cand_pos_t end) {
    std::vector<cand_pos_t> seq;

    if (end >= start) {
        seq.push_back(end - start + 1);
        for (cand_pos_t pos = 1; pos <= end - start + 1; pos++)
            seq.push_back(boustrophedon_at(start, end, pos));
    }

    return seq;
}

void W_final_pf::Sample_W(cand_pos_t start, cand_pos_t end, std::string &structure,
                          std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("W at %d and %d with W[j]=%f,%f\n", start, end, W[end], to_Energy(W[end], end));
    cand_pos_t j = end;
    cand_pos_t m = end;

    pf_t W_temp = 0;
    if (end > start) {
        for (; j > start; --j) {         // Moving through the unpaired bases in j
            if (tree.tree[j].pair < 0) { // Checking if j can be unpaired
                pf_t r = vrna_urn() * W[j];
                W_temp = W[j - 1] * scale[1];
                if (r > W_temp) { // Checking if our random sample means j is paired or unpaired
                    break;        // j is paired
                }
            } else {
                break; // j can't be unpaired so it must be paired
            }
        }
        if (j <= start + TURN) return; // No more base pairs can occur, but still successful

        pf_t r = vrna_urn() * (W[j] - W_temp);
        std::vector<cand_pos_t> is =
            boustrophedon(start, j - 1); // applies an alternating list so that the base pairing isn't biased to the right side
        cand_pos_t bous_n = is.size();
        pf_t qt = 0;
        cand_pos_t k = start;
        bool pseudoknot = false;
        if (tree.weakly_closed(1, j)) {
            for (m = 1; m < bous_n; ++m) {
                k = is[m];
                if (tree.weakly_closed(1, k - 1)) {
                    pf_t acc = (k > 1) ? W[k - 1] : 1;
                    pf_t Wkl = acc * get_energy(k, j) * exp_Extloop(k, j);
                    qt += Wkl;
                    if (qt > r) {
                        break; // k pairs with j
                    }

                    if (k == 1 || tree.weakly_closed(k, j)) {
                        Wkl = acc * get_energy_WMB(k, j) * expPS_penalty;
                        qt += Wkl;
                        if (qt > r) {
                            pseudoknot = true;
                            break; // k pairs with j as a pseudoknot
                        }
                    }
                }
            }
        }
        if (k + start > j) {
            printf("backtracking failed in ext loop at %d and %d with W[j] = %f, qt:%f < r:%f\n", start, end, W[j], qt, r);
            exit(0); /* error */
        }
        Sample_W(start, k - 1, structure, samples, tree);
        if (!pseudoknot) {
            Sample_V(k, j, structure, samples, tree);
        } else {
            Sample_WMB(k, j, structure, samples, tree);
        }
    }
}

void W_final_pf::Sample_V(cand_pos_t i, cand_pos_t j, std::string &structure,
                          std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("V at %d and %d\n", i, j);

    cand_pos_t k = i;
    cand_pos_t l = j;
    structure[i - 1] = '(';
    structure[j - 1] = ')';

    pf_t qbr = get_energy(i, j);
    pf_t V_temp = 0;

    std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
    std::pair<cand_pos_tu, cand_pos_tu> base_pair_reversed(j, i);
    ++samples[base_pair]; // Increments the base pair found in V
    ++samples[base_pair_reversed];

    pf_t r = vrna_urn() * qbr;
    pf_t qbt1 = 0;
    bool canH = !(tree.up[j - 1] < (j - i - 1));

    if (canH) V_temp = HairpinE(i, j);

    qbt1 += V_temp;
    if (qbt1 >= r) return;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1); // i+1+tree.up[i+1]?
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (k = i + 1; k <= max_k; k++) {
        if (tree.up[k - 1] >= (k - i - 1)) {
            cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
            for (l = j - 1; l >= min_l; --l) {
                if (tree.up[j - 1] >= (j - l - 1)) {
                    cand_pos_t u1 = k - i - 1;
                    cand_pos_t u2 = j - l - 1;
                    V_temp = get_energy(k, l)
                             * exp_E_IntLoop(u1, u2, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1],
                                             exp_params_);
                    V_temp *= scale[u1 + u2 + 2];
                    qbt1 += V_temp;
                    if (qbt1 >= r) break;
                }
            }
            if (qbt1 >= r) break;
        }
    }
    if (qbt1 >= r) {
        Sample_V(k, l, structure, samples, tree); // Backtrack the internal loop
        return;
    }

    V_temp = get_energy_VM(i, j); // VM includes everything since it includes the basepair (i.e. not like WM2 region), so is this fine?
    qbt1 += V_temp;
    if (qbt1 < r) {
        printf("Backtracking failed for pair (%d,%d)\n", i, j);
        exit(0);
    }

    // Must be a multiloop
    Sample_VM(i, j, structure, samples, tree);
}

void W_final_pf::Sample_VM(cand_pos_t i, cand_pos_t j, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("VM at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t qt = 0;
    if ((i + 1) + 2 * TURN + 2 >= (j - 1)) {
        printf("backtracking impossible for VM[%d, %d]\n", i, j);
        exit(0); /* error */
    }
    pf_t V_temp = 0.;
    pf_t VM_inside = get_energy_VM(i, j) / scale[2]; // If I remove scale from VM's saved values, I may save time here.
    pf_t r = vrna_urn() * VM_inside;
    bool unpaired = false;
    bool pseudoknot = false;
    for (k = i + 1; k <= j - TURN - 1; ++k) {
        V_temp = get_energy_WM(i + 1, k - 1) * get_energy_WMv(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing;
        qt += V_temp;
        if (qt > r) {
            break;
        }

        V_temp = (get_energy_WM(i + 1, k - 1) * get_energy_WMp(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        qt += V_temp;
        if (qt > r) {
            pseudoknot = true;
            break;
        }

        if (tree.up[k - 1] >= (k - (i + 1))) {
            V_temp = (expMLbase[k - i - 1] * get_energy_WMp(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
            qt += V_temp;
            if (qt > r) {
                unpaired = true;
                pseudoknot = true;
                break;
            }
        }
    }
    if (k > j - TURN) {
        printf("backtracking failed for VM at i=%d and j =%d\n", i, j);
        exit(0);
    }

    if (!unpaired) {
        Sample_WM(i + 1, k - 1, structure, samples, tree);
    }
    if (!pseudoknot) { // Case 1
        Sample_WMV(k, j - 1, structure, samples, tree);
    } else { // Case 2 or 3
        Sample_WMP(k, j - 1, structure, samples, tree);
    }
}
void W_final_pf::Sample_WM(cand_pos_t i, cand_pos_t j, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WM at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t qt = 0;
    pf_t qbt1 = 0;
    pf_t qbt2 = 0;
    bool unpaired = false;
    bool pseudoknot = false;

    pf_t V_temp = 0.;

    if (i + TURN >= j) {
        // return;
        printf("backtracking impossible for WM[%u, %u]\n", i, j);
        exit(0); /* error */
    }

    for (; j > i + TURN; --j) {
        if (tree.tree[j].pair < 0) {
            pf_t r = vrna_urn() * (get_energy_WM(i, j));

            V_temp = get_energy_WM(i, j - 1) * expMLbase[1];
            qt = V_temp;
            if (r > qt) {
                break;
            }
        } else {
            break; // j can't be unpaired so it must be paired
        }
    }

    if (i + TURN == j) {
        printf("backtracking failed for WM\n");
        exit(0); /* error */
    }

    qt = 0.;
    pf_t qm_rem = get_energy_WM(i, j) - V_temp;
    pf_t r = vrna_urn() * qm_rem;
    for (k = i; k < j - TURN; ++k) {
        qbt1 = get_energy(k, j) * exp_MLstem(k, j);
        qbt2 = get_energy_WMB(k, j) * expPSM_penalty * expb_penalty;
        bool can_pair = tree.up[k - 1] >= (k - i);
        if (can_pair) {

            V_temp = static_cast<pf_t>(expMLbase[k - i]) * qbt1;
            qt += V_temp;
            if (qt >= r) {
                unpaired = true;
                break;
            }

            V_temp = static_cast<pf_t>(expMLbase[k - i]) * qbt2;
            qt += V_temp;
            if (qt >= r) {
                unpaired = true;
                pseudoknot = true;
                break;
            }
		}

		V_temp = get_energy_WM(i, k - 1) * qbt1;
		qt += V_temp;
		if (qt >= r) break;

		V_temp = get_energy_WM(i, k - 1) * qbt2;
		qt += V_temp;
		if (qt >= r) {
			pseudoknot = true;
			break;
		}
        
    }
    if (k > j - TURN || qt < r) {
        printf("backtracking failed for WM at i=%d and j =%d with k=%d, qt=%f and r =%f and qt<r=%d\n", i, j,k,qt,r,qt<r);
        exit(0);
    }
    if (!unpaired) {
        Sample_WM(i, k - 1, structure, samples, tree);
    }
    if (!pseudoknot) {
        Sample_V(k, j, structure, samples, tree);
    } else {
        Sample_WMB(k, j, structure, samples, tree);
    }
}
void W_final_pf::Sample_WMV(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WMv at %d and %d\n", i, j);
    pf_t qt = 0;

    pf_t V_temp = 0.;

    for (; j > i + TURN; --j) {
        if (tree.tree[j].pair < 0) { // Checking if j can be unpaired
            pf_t r = vrna_urn() * get_energy_WMv(i, j);

            V_temp = get_energy_WMv(i, j - 1) * expMLbase[1];
            qt = V_temp;
            if (r > qt) {
                break;
            }
        } else {
            break; // j can't be unpaired so it must be paired
        }
    }

    if (i + TURN == j) {
        printf("backtracking failed for WMV\n");
        exit(0); /* error */
    }

    Sample_V(i, j, structure, samples, tree);
}

void W_final_pf::Sample_WMP(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WMp at %d and %d\n", i, j);
    pf_t qt = 0;

    pf_t V_temp = 0.;

    for (; j > i + TURN; --j) {
        if (tree.tree[j].pair < 0) { // Checking if j can be unpaired
            pf_t r = vrna_urn() * get_energy_WMp(i, j);

            V_temp = get_energy_WMp(i, j - 1) * expMLbase[1];
            qt = V_temp;
            if (r > qt) {
                break;
            }
        } else {
            break; // j can't be unpaired so it must be paired
        }
    }

    if (i + TURN == j) { // I'm kinda assuming something like ([..)] for this
        printf("backtracking failed for WMP\n");
        exit(0); /* error */
    }

    Sample_WMB(i, j, structure, samples, tree);
}

void W_final_pf::Sample_WMB(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WMB at %d and %d\n", i, j);
    cand_pos_t l = j;
    pf_t qt = 0;
    cand_pos_t Bp_lj = 0;
    cand_pos_t bp_j = 0;

    pf_t V_temp = 0.;

    pf_t r = vrna_urn() * get_energy_WMB(i, j);

    if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i) {
        bp_j = tree.tree[j].pair;
        for (l = (bp_j + 1); (l < j); ++l) {
            // if(tree.tree[l].pair>0) continue;
            Bp_lj = tree.Bp(l, j);
            if (Bp_lj >= 0 && Bp_lj < n) {
                V_temp = get_BE(bp_j, j, tree.tree[Bp_lj].pair, Bp_lj, tree) * get_energy_WMBP(i, l) * get_energy_WI(l + 1, Bp_lj - 1) * expPB_penalty;
                qt += V_temp;
                if (qt > r) {
                    break;
                }
            }
        }
    }
    if (qt >= r) { // I could put this in the for loop then just do sample_WMBP if it doesn't sample in there
        Sample_BE(bp_j, j, tree.tree[Bp_lj].pair, Bp_lj, structure, samples, tree);
        Sample_WMBP(i, l, structure, samples, tree);
        Sample_WI(l + 1, Bp_lj - 1, structure, samples, tree);
		return;
    }
    V_temp = get_energy_WMBP(i, j);
    qt += V_temp;
    if (qt <= r) {
        printf("backtracking failed for WMB\n");
        exit(0); /* error */
    }
    Sample_WMBP(i, j, structure, samples, tree);
}

void W_final_pf::Sample_WI(cand_pos_t i, cand_pos_t j, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WI at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t qt = 0, qbt1 = 0, qbt2 = 0;
    bool pseudoknot = false;

    pf_t V_temp = 0;
    if (j > i) {
        for (; j > i + TURN; --j) {
            if (tree.tree[j].pair < 0) { // Checking if j can be unpaired
                pf_t r = vrna_urn() * (get_energy_WI(i, j));

                V_temp = get_energy_WI(i, j - 1) * expPUP_pen[1];
                qt = V_temp;
                if (r > qt) {
                    break;
                }
            } else {
                break; // j can't be unpaired so it must be paired
            }
        }
        if (j <= i + TURN) return; // No more base pairs can occur, but still successful

        qt = 0;
        pf_t qm_rem = get_energy_WI(i, j) - V_temp;

        pf_t r = vrna_urn() * qm_rem;

        for (k = i; k <= j - TURN - 1; k++) {
            qbt1 = get_energy(k, j) * expPPS_penalty;
            qbt2 = get_energy_WMB(k, j) * expPSP_penalty * expPPS_penalty;

            V_temp = qbt1 * get_energy_WI(i, k - 1);
            qt += V_temp;
            if (qt >= r) break;

            V_temp = qbt2 * get_energy_WI(i, k - 1);
            qt += V_temp;
            if (qt >= r) {
                pseudoknot = true;
                break;
            }
        }

        Sample_WI(i, k - 1, structure, samples, tree);
        if (!pseudoknot) {
            Sample_V(k, j, structure, samples, tree);
        } else {
            Sample_WMB(k, j, structure, samples, tree);
        }
    }
}

void W_final_pf::Sample_WIP(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WIP at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t fbd = 0;
    pf_t qt = 0, qbt1 = 0, qbt2 = 0;
    bool unpaired = false;
    bool pseudoknot = false;

    pf_t V_temp = 0;
    if (j <= i) return;
    for (; j > i + TURN; --j) {
        if (tree.tree[j].pair < 0) { // Checking if j can be unpaired
            pf_t r = vrna_urn() * (get_energy_WIP(i, j) - fbd);

            V_temp = get_energy_WIP(i, j - 1) * expcp_pen[1];
            qt = V_temp;
            if (r > qt) {
                break;
            }
        } else {
            break; // j can't be unpaired so it must be paired
        }
    }
    if (i + TURN == j) {
        printf("backtracking failed for WIP\n");
        exit(0); /* error */
    }

    qt = 0;
    pf_t qm_rem = get_energy_WIP(i, j) - V_temp;

    pf_t r = vrna_urn() * (qm_rem - fbd);

    for (k = i; k < j - TURN; ++k) {
        qbt1 = get_energy(k, j) * expbp_penalty;
        qbt2 = get_energy_WMB(k, j) * expbp_penalty * expPSM_penalty;

        bool can_pair = tree.up[k - 1] >= (k - i);
        if (can_pair) {
            V_temp = qbt1 * expcp_pen[k - i];
            qt += V_temp;
            if (qt >= r) {
                unpaired = true;
                break;
            }

            V_temp = qbt2 * expcp_pen[k - i];
            qt += V_temp;
            if (qt >= r) {
                unpaired = true;
                pseudoknot = true;
                break;
            }
        }

        if (k >= i) { // It had > but for me, it should be >=?
            V_temp = qbt1 * get_energy_WM(i, k - 1);
            qt += V_temp;
            if (qt >= r) break;

            V_temp = qbt2 * get_energy_WM(i, k - 1);
            qt += V_temp;
            if (qt >= r) {
                pseudoknot = true;
                break;
            }
        }
    }

    if (k + TURN >= j) {
        printf("backtracking failed for WIP right base pair with k=%d and j =%d, and qt=%f with r-%f\n", k, j, qt, r);
        exit(0);
    }
    if (!unpaired) {
        Sample_WIP(i, k - 1, structure, samples, tree);
    }
    if (!pseudoknot) {
        Sample_V(k, j, structure, samples, tree);
    } else {
        Sample_WMB(k, j, structure, samples, tree);
    }
}

void W_final_pf::Sample_WMBW(cand_pos_t i, cand_pos_t j, std::string &structure,
                             std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WMBW at %d and %d\n", i, j);
    cand_pos_t l = j;
    pf_t fbd = 0;
    pf_t qt = 0;

    pf_t V_temp = 0;
    pf_t r = vrna_urn() * (get_energy_WMBW(i, j) - fbd);
    if (tree.tree[j].pair < j) {
        for (l = i + 1; l < j; l++) {
            if (tree.tree[l].pair < 0 && tree.tree[l].parent->index > -1 && tree.tree[j].parent->index > -1
                && tree.tree[j].parent->index == tree.tree[l].parent->index) {
                V_temp = get_energy_WMBP(i, l) * get_energy_WI(l + 1, j);
                qt += V_temp;
                if (qt >= r) {
                    break;
                }
            }
        }
    }
    if (l >= j) {
        printf("Backtracking failed in WMBW for pair (%d,%d) with qt=%f < r=%f and l = %d\n", i, j, qt, r, l);
        exit(0);
    }
    Sample_WMBP(i, l, structure, samples, tree);
    Sample_WI(l + 1, j, structure, samples, tree);
	return;
}

void W_final_pf::Sample_WMBP(cand_pos_t i, cand_pos_t j, std::string &structure,
                             std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("WMBP at %d and %d\n", i, j);
    cand_pos_t l = j;
    pf_t qt = 0;
    cand_pos_t bp_il = 0;
    cand_pos_t Bp_lj = 0;
    cand_pos_t B_lj = 0;
    cand_pos_t b_ij = tree.b(i, j);
    bool case1 = false, case2 = false, case4 = false;

    pf_t V_temp = 0;
    pf_t r = vrna_urn() * get_energy_WMBP(i, j);

    if (tree.tree[j].pair < 0) {
        for (l = i + 1; l < j - TURN; ++l) {
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j, tree);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                bp_il = tree.bp(i, l);
                Bp_lj = tree.Bp(l, j);
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) {
                    B_lj = tree.B(l, j);
                    if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l + TURN <= j) {
                        V_temp = get_BE(tree.tree[B_lj].pair, B_lj, tree.tree[Bp_lj].pair, Bp_lj, tree) * get_energy_WMBP(i, l - 1)
                                 * get_energy_VP(l, j) * pow(expPB_penalty, 2);
                        qt += V_temp;
                        if (qt >= r) {
                            case1 = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (case1) {
        Sample_BE(tree.tree[B_lj].pair, B_lj, tree.tree[Bp_lj].pair, Bp_lj, structure, samples, tree);
        Sample_WMBP(i, l - 1, structure, samples, tree);
        Sample_VP(l, j, structure, samples, tree);
        return;
    }

    if (tree.tree[j].pair < 0) {
        for (l = i + 1; l < j - TURN; l++) {
            bp_il = tree.bp(i, l);
            Bp_lj = tree.Bp(l, j);
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(l, j, tree);
            if ((b_ij > 0 && l < b_ij) || (b_ij < 0 && ext_case == 0)) {
                if (bp_il >= 0 && l > bp_il && Bp_lj > 0 && l < Bp_lj) {
                    B_lj = tree.B(l, j);
                    if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l + TURN <= j) {
                        V_temp = get_BE(tree.tree[B_lj].pair, B_lj, tree.tree[Bp_lj].pair, Bp_lj, tree) * get_energy_WMBW(i, l - 1)
                                 * get_energy_VP(l, j) * pow(expPB_penalty, 2);
                        qt += V_temp;
                        if (qt >= r) {
                            case2 = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (case2) {
        Sample_BE(tree.tree[B_lj].pair, B_lj, tree.tree[Bp_lj].pair, Bp_lj, structure, samples, tree);
        Sample_WMBW(i, l - 1, structure, samples, tree);
        Sample_VP(l, j, structure, samples, tree);
        return;
    }

    V_temp = get_energy_VP(i, j) * expPB_penalty;
    qt += V_temp;
    if (qt >= r) {
        Sample_VP(i, j, structure, samples, tree);
        return;
    }

    if (tree.tree[j].pair < 0 && tree.tree[i].pair >= 0) {
        for (l = i + 1; l < j; l++) {
            bp_il = tree.bp(i, l);
            if (bp_il >= 0 && bp_il < n && l + TURN <= j) {
                if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l + TURN <= j) {
                    V_temp = get_BE(i, tree.tree[i].pair, bp_il, tree.tree[bp_il].pair, tree) * get_energy_WI(bp_il + 1, l - 1) * get_energy_VP(l, j)
                             * pow(expPB_penalty, 2);
                    qt += V_temp;
                    if (qt >= r) {
                        case4 = true;
                        break;
                    }
                }
            }
        }
    }
    if (case4) {
        Sample_BE(i, tree.tree[i].pair, bp_il, tree.tree[bp_il].pair, structure, samples, tree);
        Sample_WI(bp_il + 1, l - 1, structure, samples, tree);
        Sample_VP(l, j, structure, samples, tree);
		return;
    } else {
        printf("backtracking failed for WMBP\n");
        exit(0);
    }
}

void W_final_pf::Sample_VP(cand_pos_t i, cand_pos_t j, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("VP at %d and %d\n", i, j);
    cand_pos_t k, l;
    pf_t qt = 0;
    structure[i - 1] = '[';
    structure[j - 1] = ']';

    pf_t V_temp = 0;
    pf_t r = vrna_urn() * get_energy_VP(i, j);

    std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
    std::pair<cand_pos_tu, cand_pos_tu> base_pair_reversed(j, i);
    ++samples[base_pair]; // Increments the base pair found in VP
    ++samples[base_pair_reversed];

    cand_pos_t Bp_ij = tree.Bp(i, j);
    cand_pos_t B_ij = tree.B(i, j);
    cand_pos_t b_ij = tree.b(i, j);
    cand_pos_t bp_ij = tree.bp(i, j);

    if ((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) < (tree.tree[i].parent->index) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        V_temp = (get_energy_WI(i + 1, Bp_ij - 1) * get_energy_WI(B_ij + 1, j - 1));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt >= r) {
            Sample_WI(i + 1, Bp_ij - 1, structure, samples, tree);
            Sample_WI(B_ij + 1, j - 1, structure, samples, tree);
            return;
        }
    }
    if ((tree.tree[i].parent->index) < (tree.tree[j].parent->index) && (tree.tree[j].parent->index) > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        V_temp = (get_energy_WI(i + 1, b_ij - 1) * get_energy_WI(bp_ij + 1, j - 1));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt >= r) {
            Sample_WI(i + 1, b_ij - 1, structure, samples, tree);
            Sample_WI(bp_ij + 1, j - 1, structure, samples, tree);
            return;
        }
    }
    if ((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        V_temp = (get_energy_WI(i + 1, Bp_ij - 1) * get_energy_WI(B_ij + 1, b_ij - 1) * get_energy_WI(bp_ij + 1, j - 1));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt >= r) {
            Sample_WI(i + 1, Bp_ij - 1, structure, samples, tree);
            Sample_WI(B_ij + 1, b_ij - 1, structure, samples, tree);
            Sample_WI(bp_ij + 1, j - 1, structure, samples, tree);
            return;
        }
    }
    pair_type ptype_closingip1jm1 = pair[S_[i + 1]][S_[j - 1]];
    if ((tree.tree[i + 1].pair) < -1 && (tree.tree[j - 1].pair) < -1 && ptype_closingip1jm1 > 0) {
        V_temp = (get_e_stP(i, j) * get_energy_VP(i + 1, j - 1));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt >= r) {
            Sample_VP(i + 1, j - 1, structure, samples, tree);
            return;
        }
    }

    cand_pos_t min_borders = std::min((cand_pos_tu)Bp_ij, (cand_pos_tu)b_ij);
    cand_pos_t edge_i = std::min(i + MAXLOOP + 1, j - TURN - 1);
    min_borders = std::min(min_borders, edge_i);
    for (k = i + 1; k < min_borders; ++k) {
        if (tree.tree[k].pair < -1 && (tree.up[(k)-1] >= ((k) - (i)-1))) {
            cand_pos_t max_borders = std::max(bp_ij, B_ij) + 1;
            cand_pos_t edge_j = k + j - i - MAXLOOP - 2;
            max_borders = std::max(max_borders, edge_j);
            for (l = j - 1; l > max_borders; --l) {
                pair_type ptype_closingkj = pair[S_[k]][S_[l]];
                if (k == i + 1 && l == j - 1) continue; // I have to add or else it will add a stP version and an eintP version to the sum
                if (tree.tree[l].pair < -1 && ptype_closingkj > 0 && (tree.up[(j)-1] >= ((j) - (l)-1))) {
                    cand_pos_t u1 = k - i - 1;
                    cand_pos_t u2 = j - l - 1;
                    V_temp = (get_e_intP(i, k, l, j) * get_energy_VP(k, l));
                    V_temp *= scale[u1 + u2 + 2];
                    qt += V_temp;
                    if (qt >= r) {
                        break;
                    }
                }
            }
            if (qt >= r) {
                break;
            }
        }
    }
    if (k < min_borders) {
        Sample_VP(k, l, structure, samples, tree);
        return;
    }

    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree.b(i, j), (cand_pos_tu)tree.Bp(i, j));
    cand_pos_t max_i_bp = std::max(tree.B(i, j), tree.bp(i, j));
    for (k = i + 1; k < min_Bp_j; ++k) {
        V_temp = (get_energy_WIP(i + 1, k - 1) * get_energy_VP(k, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt > r) {
            break;
        }
    }
    if (k < min_Bp_j) {
        Sample_WIP(i + 1, k - 1, structure, samples, tree);
        Sample_VP(k, j - 1, structure, samples, tree);
		return;
    }

    for (k = max_i_bp + 1; k < j; ++k) {
        V_temp = (get_energy_VP(i + 1, k) * get_energy_WIP(k + 1, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt > r) {
            break;
        }
    }
    if (k < j) {
        Sample_VP(i + 1, k, structure, samples, tree);
        Sample_WIP(k + 1, j - 1, structure, samples, tree);
        return;
    }

    for (k = i + 1; k < min_Bp_j; ++k) {
        V_temp = (get_energy_WIP(i + 1, k - 1) * get_energy_VPR(k, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt > r) {
            break;
        }
    }
    if (k < min_Bp_j) {
        Sample_WIP(i + 1, k - 1, structure, samples, tree);
        Sample_VPR(k, j - 1, structure, samples, tree);
        return;
    }

    for (k = max_i_bp + 1; k < j; ++k) {
        V_temp = (get_energy_VPL(i + 1, k) * get_energy_WIP(k + 1, j - 1) * expap_penalty * pow(expbp_penalty, 2));
        V_temp *= scale[2];
        qt += V_temp;
        if (qt > r) {
            break;
        }
    }
    if (k < j) {
        Sample_VPL(i + 1, k, structure, samples, tree);
        Sample_WIP(k + 1, j - 1, structure, samples, tree);
        return;
    }
}

void W_final_pf::Sample_VPL(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("VPL at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t fbd = 0;
    pf_t qt = 0;

    pf_t V_temp = 0;
    pf_t r = vrna_urn() * (get_energy_VPL(i, j) - fbd);
    cand_pos_t min_Bp_j = std::min((cand_pos_tu)tree.b(i, j), (cand_pos_tu)tree.Bp(i, j));
    for (k = i + 1; k < min_Bp_j; ++k) {
        bool can_pair = tree.up[k - 1] >= (k - i);
        if (can_pair) {
            V_temp = (expcp_pen[k - i] * get_energy_VP(k, j));
            qt += V_temp;
            if (qt >= r) {
                break;
            }
        }
    }
    if (k < min_Bp_j) {
        Sample_VP(k, j, structure, samples, tree);
    } else {
        printf("Backtracking error in VPL\n");
        exit(0);
    }
}

void W_final_pf::Sample_VPR(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    if (debug) printf("VPR at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t fbd = 0;
    pf_t qt = 0;
    bool unpaired = false;

    pf_t V_temp = 0;
    pf_t r = vrna_urn() * (get_energy_VPR(i, j) - fbd);

    cand_pos_t max_i_bp = std::max(tree.B(i, j), tree.bp(i, j));
    for (k = max_i_bp + 1; k < j; ++k) {
        bool can_pair = tree.up[j - 1] >= (j - k);
        V_temp = (get_energy_VP(i, k) * get_energy_WIP(k + 1, j));
        qt += V_temp;
        if (qt >= r) {
            break;
        }
        if (can_pair) {
            V_temp = (get_energy_VP(i, k) * expcp_pen[k - i]);
            qt += V_temp;
            if (qt >= r) {
                unpaired = true;
                break;
            }
        }
    }
    if (k == j) {
        printf("Backtracking error in VPR\n");
        exit(0);
    }

    if (!unpaired) {
        Sample_WIP(k + 1, j, structure, samples, tree);
    }
    Sample_VP(i, k, structure, samples, tree);
}

void W_final_pf::Sample_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {
    
	if (debug) printf("BE at %d and %d, and %d and %d\n", i, j, ip, jp);

    if (!(i >= 1 && i <= ip && ip < jp && jp <= j && j <= n && tree.tree[i].pair > 0 && tree.tree[j].pair > 0 && tree.tree[ip].pair > 0
          && tree.tree[jp].pair > 0)) { // impossible cases
        printf("Backtracking failed in BE: impossible case -- %d and %d, and %d and %d\n", i, j, ip, jp);
        exit(0);
    }

    if (tree.tree[i].pair != j || tree.tree[ip].pair != jp) {
        printf("Backtracking failed in BE: base case: i.j and ip.jp must be in G\n");
        exit(0);
    }

	cand_pos_t l = j;
    cand_pos_t lp = j;
    pf_t qt = 0;
    bool unpaired_left = false;
    bool unpaired_right = false;
    structure[i - 1] = '(';
    structure[j - 1] = ')';
    std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
    std::pair<cand_pos_tu, cand_pos_tu> base_pair_reversed(j, i);
    ++samples[base_pair]; // Increments the base pair found in BE
    ++samples[base_pair_reversed];

    pf_t V_temp = 0;
    pf_t r = vrna_urn() * get_BE(i, j, ip, jp, tree);

    if (i == ip && j == jp && i < j) {
        return;
    }

    if (tree.tree[i + 1].pair == j - 1) {
        Sample_BE(i + 1, j - 1, ip, jp, structure, samples, tree);
        return;
    }

    for (l = i + 1; l <= ip; l++) {
        if (tree.tree[l].pair >= -1 && jp <= tree.tree[l].pair && tree.tree[l].pair < j) {
            lp = tree.tree[l].pair;

            bool empty_region_il = (tree.up[(l)-1] >= l - i - 1);       // empty between i+1 and l-1
            bool empty_region_lpj = (tree.up[(j)-1] >= j - lp - 1);     // empty between lp+1 and j-1
            bool weakly_closed_il = tree.weakly_closed(i + 1, l - 1);   // weakly closed between i+1 and l-1
            bool weakly_closed_lpj = tree.weakly_closed(lp + 1, j - 1); // weakly closed between lp+1 and j-1
			pf_t expbp2 = pow(expbp_penalty, 2);

            if (empty_region_il && empty_region_lpj) {
                cand_pos_t u1 = l - i - 1;
                cand_pos_t u2 = j - lp - 1;
                V_temp = get_e_intP(i, l, lp, j) * get_BE(l, lp, ip, jp, tree);
                V_temp *= scale[u1 + u2 + 2];
                qt += V_temp; // Added to e_intP that l != i+1 and lp != j-1 at the same time
                if (qt >= r) {
                    unpaired_left = true;
                    unpaired_right = true;
                    break;
                }
            }
            if (weakly_closed_il && weakly_closed_lpj) {
                V_temp = get_energy_WIP(i + 1, l - 1) * get_BE(l, lp, ip, jp, tree) * get_energy_WIP(lp + 1, j - 1) * expap_penalty * expbp2;
                V_temp *= scale[2];
                qt += V_temp;
                if (qt >= r) {
                    break;
                }
            }
            if (weakly_closed_il && empty_region_lpj) {
                V_temp = get_energy_WIP(i + 1, l - 1) * get_BE(l, lp, ip, jp, tree) * expcp_pen[j - lp - 1] * expap_penalty * expbp2;
                V_temp *= scale[2];
                qt += V_temp;
                if (qt >= r) {
                    unpaired_right = true;
                    break;
                }
            }
            if (empty_region_il && weakly_closed_lpj) {
                V_temp = expcp_pen[l - i - 1] * get_BE(l, lp, ip, jp, tree) * get_energy_WIP(lp + 1, j - 1) * expap_penalty * expbp2;
                V_temp *= scale[2];
                qt += V_temp;
                if (qt >= r) {
                    unpaired_left = true;
                    break;
                }
            }
        }
    }

    if (!unpaired_left) {
        Sample_WIP(i + 1, l - 1, structure, samples, tree);
    }
    if (!unpaired_right) {
        Sample_WIP(lp + 1, j - 1, structure, samples, tree);
    }
    Sample_BE(l, lp, ip, jp, structure, samples, tree);
}
char W_final_pf::bpp_symbol(cand_pos_t i, cand_pos_t j, pf_t *P, sparse_tree &tree) {
    if (i > j) {
        cand_pos_t temp = i;
        i = j;
        j = temp;
    }
    bool weakly_closed_ij = tree.weakly_closed(i, j);
    if (P[0] > 0.667) return '.';
    if (P[0] > (P[1] + P[2])) return ',';
    if (weakly_closed_ij) {
        if (P[1] > 0.667) return '(';
        if (P[2] > 0.667) return ')';
        if ((P[1] + P[2]) > P[0]) {
            if ((P[1] / (P[1] + P[2])) > 0.667) return '{';

            if ((P[2] / (P[1] + P[2])) > 0.667)
                return '}';
            else
                return '|';
        }
    } else {
        if (P[1] > 0.667) return '[';
        if (P[2] > 0.667) return ']';
        if ((P[1] + P[2]) > P[0]) {
            if ((P[1] / (P[1] + P[2])) > 0.667) return '/';

            if ((P[2] / (P[1] + P[2])) > 0.667)
                return '\\';
            else
                return '|';
        }
    }
    return ':';
}
void W_final_pf::pairing_tendency(std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples, sparse_tree &tree) {

    for (cand_pos_t j = 1; j <= n; j++) {
        cand_pos_t best_i = j;
        pf_t best_frequency = 0;
        pf_t P[3] = {1, 0, 0};
        for (cand_pos_t i = 1; i < j; i++) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            pf_t frequency_ij = (pf_t)samples[base_pair] / num_samples;
            P[2] += frequency_ij;
            P[0] -= frequency_ij;
            if (frequency_ij > best_frequency) {
                best_frequency = frequency_ij;
                best_i = i;
            }
        }
        for (cand_pos_t i = j + 1; i <= n; i++) {
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(j, i);
            pf_t frequency_ji = (pf_t)samples[base_pair] / num_samples;
            P[1] += frequency_ji;
            P[0] -= frequency_ji;
            if (frequency_ji > best_frequency) {
                best_frequency = frequency_ji;
                best_i = i;
            }
        }
        structure[j - 1] = bpp_symbol(best_i, j, P, tree);
    }
}