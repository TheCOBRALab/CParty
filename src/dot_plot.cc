#include "dot_plot.hh"
#include <iostream>
#include <time.h>

void create_PS_header(std::ofstream &out){
    out << "%%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
    out << "%%%%Creator: COBRAlab" << std::endl;
    out << "%%%%CreationDate: " << std::endl;
    out << "%%%%Title: Dot.ps" << std::endl;
    out << "%%%%BoundingBox: 66 211 518 662" << std::endl;
    out << "%%%%DocumentFonts: Helvetica" << std::endl;
    out << "%%%%Pages: 1" << std::endl;
    out << "%%%%EndComments" << std::endl;
    out << std::endl;

    out << "/DPdict 100 dict def" << std::endl;
    out << std::endl;

    out << "DPdict begin" << std::endl;
    out << std::endl;
    out << "%%%%BeginProlog" << std::endl;
    out << std::endl;

    out << PS_dot_plot_base << std::endl;
    out << PS_dot_plot_sd << std::endl;
    out << PS_dot_plot_ud << std::endl;
    out << PS_dot_plot_sc_motifs << std::endl;
    out << PS_dot_plot_linear_data << std::endl;

    out << "%%%%EndProlog" << std::endl;
    out << std::endl;

}
void create_PS_title(std::ofstream &out){
    out << "/DPtitle {\n  (%s)\n} def" << std::endl;
    out << std::endl;
}

void create_PS_sequence(std::ofstream &out, std::string &seq){
    out << "/sequence { (\\" <<std::endl;
    out << seq << std::endl;
    out << ") } def" <<std::endl;
    out << "/len { sequence length } bind def" <<std::endl;
    out << std::endl;
}

void create_PS_footer(std::ofstream &out){
    out << "showpage" << std::endl;
    out << "end" << std::endl;
    out << "%%%%EOF" << std::endl;
}

void create_PS_data(std::ofstream &out,std::unordered_map< std::pair<cand_pos_t,cand_pos_t>,cand_pos_t, SzudzikHash > &samples, int num_samples, cand_pos_t n){
    int min_samples = .1*num_samples;
    for(cand_pos_t i = 1; i <=n; ++i){
        for(cand_pos_t j = 1; j <=n; ++j){
            std::pair <cand_pos_tu,cand_pos_tu> base_pair (i,j);
            if(samples[base_pair]>min_samples){
                pf_t prob = (pf_t) samples[base_pair]/num_samples;
                out << i << " " << j << " " << prob << " ubox" << std::endl;
            }

        }
    }
    for(cand_pos_t i = 1; i <=n; ++i){
        for(cand_pos_t j = 1; j <=n; ++j){
            std::pair <cand_pos_tu,cand_pos_tu> base_pair (i,j);
            if(samples[base_pair]>min_samples){
                pf_t prob = (pf_t) samples[base_pair]/num_samples;
                out << i << " " << j << " " << prob << " lbox" << std::endl;
            }

        }
    }
}

void create_dot_plot(std::string &seq, std::unordered_map< std::pair<cand_pos_t,cand_pos_t>,cand_pos_t, SzudzikHash > &samples, int num_samples){

    std::ofstream out("Dot.ps");
    cand_pos_t n = seq.length();

    create_PS_header(out);
    create_PS_title(out);
    create_PS_sequence(out,seq);

    out << "72 216 translate\n72 6 mul len 1 add div dup scale" << std::endl;
    out << "/Helvetica findfont 0.95 scalefont setfont" << std::endl;
    out << std::endl;

    out << "drawseq" <<std::endl;
    out << "%%data starts here" << std::endl;
    out << std::endl;

    out << "%%draw the grid\ndrawgrid" << std::endl;
    out <<std::endl;
    out << "%%start of base pair probability data" << std::endl;
    
    create_PS_data(out,samples,num_samples,n);
    create_PS_footer(out);
    out.close();

}