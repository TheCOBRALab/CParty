#ifndef H_dotplot_H
#define H_dotplot_H

#include "part_func.hh"
#include <sparse_tree.hh>
#include "base_types.hh"
#include <math.h>
#include <fstream>

static const unsigned char PS_dot_plot_base[] = {
    #include "postscript/dot_plot_base.hex"
};

static const unsigned char PS_dot_plot_sd[] = {
    #include "postscript/dot_plot_sd.hex"
};

static const unsigned char PS_dot_plot_ud[] = {
    #include "postscript/dot_plot_ud.hex"
};

static const unsigned char PS_dot_plot_sc_motifs[] = {
    #include "postscript/dot_plot_sc_motifs.hex"
};

static const unsigned char PS_dot_plot_linear_data[] = {
    #include "postscript/dot_plot_linear_data.hex"
};

void create_PS_header(std::ofstream &out);

void create_PS_title(std::ofstream &out);

void create_PS_sequence(std::ofstream &out, std::string &seq);

void create_PS_footer(std::ofstream &out);

void create_PS_data(std::ofstream &out,std::unordered_map< std::pair<cand_pos_t,cand_pos_t>,cand_pos_t, SzudzikHash > &samples,std::vector<Node> tree,std::string MFE_structure, int num_samples, cand_pos_t n);

void create_dot_plot(std::string &seq,std::vector<Node> tree,std::string &MFE_structure, std::unordered_map< std::pair<cand_pos_t,cand_pos_t>,cand_pos_t, SzudzikHash > &samples, int num_samples);


#endif