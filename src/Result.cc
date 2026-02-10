#include "Result.hh"

// constructor
Result::Result(std::string sequence, std::string restricted, double restricted_energy, std::string final_structure, double final_energy,
               std::string final_structure_pf, pf_t pf_energy, std::string MEA_structure, pf_t MEA, std::string centroid_structure, pf_t distance,std::string fatgraph,pf_t fatgraph_frquency, pf_t frequency, pf_t diversity) {
    this->sequence = sequence;
    this->restricted = restricted;
    this->restricted_energy = restricted_energy;
    this->final_structure = final_structure;
    this->final_energy = final_energy;
    this->final_structure_pf = final_structure_pf;
    this->pf_energy = pf_energy;
    this->MEA_structure = MEA_structure;
    this->MEA = MEA;
    this->centroid_structure = centroid_structure;
    this->distance = distance;
    this->fatgraph = fatgraph;
    this->fatgraph_frequency = fatgraph_frquency;
    this->frequency = frequency;
    this->diversity = diversity;
}

// destructor
Result::~Result() {}

std::string Result::get_sequence() { return this->sequence; }
std::string Result::get_restricted() { return this->restricted; }
double Result::get_restricted_energy() { return this->restricted_energy; }
std::string Result::get_final_structure() { return this->final_structure; }
double Result::get_final_energy() { return this->final_energy; }
std::string Result::get_final_structure_pf() { return this->final_structure_pf; }
pf_t Result::get_pf_energy() { return this->pf_energy; }
std::string Result::get_MEA_structure() { return this->MEA_structure; }
pf_t Result::get_MEA() { return this->MEA; }
std::string Result::get_centroid_structure() { return this->centroid_structure; }
pf_t Result::get_distance() { return this->distance; }
std::string Result::get_fatgraph() { return this->fatgraph; }
pf_t Result::get_fatgraph_frequency() { return this->fatgraph_frequency; }
pf_t Result::get_frequency() { return this->frequency; }
pf_t Result::get_diversity() { return this->diversity; }
