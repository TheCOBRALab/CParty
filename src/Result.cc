#include "Result.hh"

//constructor
Result::Result(std::string sequence,std::string restricted,double restricted_energy, std::string final_structure, double final_energy,std::string final_structure_pf, pf_t pf_energy){
    this->sequence = sequence;
    this->restricted = restricted;
    this->restricted_energy = restricted_energy;
    this->final_structure = final_structure;
    this->final_energy = final_energy;
    this->final_structure_pf = final_structure_pf;
    this->pf_energy = pf_energy;

}

//destructor
Result::~Result(){
   
}

std::string Result::get_sequence(){
    return this->sequence;
}
std::string Result::get_restricted(){
    return this->restricted;
}
std::string Result::get_final_structure(){
    return this->final_structure;
}
double Result::get_final_energy(){
    return this->final_energy;
}
std::string Result::get_final_structure_pf(){
    return this->final_structure_pf;
}

pf_t Result::get_pf_energy(){
    return this->pf_energy;
}

double Result::get_restricted_energy(){
    return this->restricted_energy;
}
