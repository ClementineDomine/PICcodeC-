/**
Written by: Daniel Duque and Clementine Domine
Last modified on 28 fev 2020

Main file for plasma simulation.
This file is used to modify the objectives of the simulation.

Currently generates the files to investigate what values of Nr give us a tolerable accuracy and convergence rate

/*
Main file for plasma simulation.
This file is used to modify the objectives of the simulation.

Currently generates the files to investigate what values of Nr give us a tolerable accuracy and convergence rate
*/
#include"PenningTrap.hpp"
#include"Plasma.hpp"
#include"Constants.hpp"
#include<vector>
#include<fstream>
#include<iomanip>
#include<iostream>

//This attempts to study convergence rate for Nr.
int main(void)
{
  int Nz{ 400 };
    int Nr{ 150 };
    double trapRadius{ 0.01488 };
    double aLength{ 0.01322 };
    double aGap{ 0.0005 };
    std::vector<Electrode> trapElectrodes;
    trapElectrodes.push_back(Electrode(aLength, -50));
    trapElectrodes.push_back(Electrode(aLength, 0));
    trapElectrodes.push_back(Electrode(aLength, 0));
    trapElectrodes.push_back(Electrode(aLength, -50));
    std::vector<double> trapGaps(3, aGap);
    PenningTrap theTrap(trapRadius, trapElectrodes, trapGaps, Nz, Nr);
    Plasma electronPlasma(theTrap, "electrons", massE, -ePos);
    electronPlasma.loadProfile(-1e-12,3.5,3.5,100);
    return 0;
}
   
    
   //   electrons.extractselfpotential11 ("./Semester 2/Potentialenergy/correctionpotential11.csv", 0.02669);
   //   electrons.extractselfpotential12 ("./Semester 2/Potentialenergy/correctionpotential12.csv", 0.02669);
  //   electrons.loadSingleRing( -100 * ePos, 0, 0.02669, 0);
    
//    electrons.loadOneDUniform(50, -100 * ePos, 0.01, 1);
//    electrons.gettotalpotentialenergy(0,0.02669);
  //   electrons.extracttotalpotentialenergy("./Semester 2//Potentialenergy/potentialenergycorrection.csv", 0,0.02669);
    // double currentTime{ 0 };
   // theTrap.saveStates(currentTime);
  //  theTrap.extractTrapParameters("./Semester 2/Potentialenergy/ TrapParameters.csv");
  //  theTrap.extractPlasmasHistories("./Semester 2/Potentialenergy/.");
   // electrons.extractPlasmaParameters("./Semester 2/Potentialenergy/ ElectronsParameters.csv");
  //   electrons.extractSelfPotential("./Semester 2/Potentialenergy/Self.csv");//Store the potential field in a document called filename
