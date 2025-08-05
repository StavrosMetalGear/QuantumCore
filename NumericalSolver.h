#pragma once
#include <vector>
#include <string>

class NumericalSolver {
public:
    static void solveSchrodingerFDM(
        double mass,
        double xMin,
        double xMax,
        int numPoints,
        const std::vector<double>& potential,
        int numEigenstates,
        const std::string& outputFilename
    );
};

