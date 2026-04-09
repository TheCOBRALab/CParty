#ifndef SHAPE
#define SHAPE

#include <vector>
#include <string>
#include "base_types.hh"

// Consider floats instead of doubles
class SHAPEData {

    private:
        std::vector<double> calculated;  // pseudo-energy calculated SHAPE values, indexed by nucleotide position
        std::vector<double> expcalculated; // Boltzmann variation of the calculated values
        cand_pos_t n;
        double  slope;         // m parameter for linear conversion to energy bonus
        double  intercept;     // b parameter

        double calculate(double reactivity);
        bool exists(const std::string &filename);
    public:
        SHAPEData(const std::string &filename, cand_pos_t n, double slope = 1.8, double intercept = -.6);
        double get_calculated(cand_pos_t i);
        double get_expcalculated(cand_pos_t i);
        void rescale_calculate(double kT, double TT, int smooth);
};


#endif