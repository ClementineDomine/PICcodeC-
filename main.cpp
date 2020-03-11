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
  
// for (int voltage=-15; voltage> -40; voltage -= 1)
    for (int i =0;i<3;i++)
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
        electronPlasma.loadProfile(1000,-1e-13,3.5,3.5,1e4);
   }
//    double Endoftimestep=10;
//    for(int deltat = 0; deltat < Endoftimestep+1; ++deltat)
//          {
//              theTrap.setPotential(4, ((-50-(deltat/Endoftimestep)*voltage)));
//              theTrap.movePlasmas(10e-12);
//          }
          //    for(int deltat = 0; deltat < 300; ++deltat)
//          {
//              theTrap.movePlasmas(10e-12);
//          }
//    double Size= electronPlasma.extractsizeringvector();
//    std::ofstream myfile;
//    myfile.open ("numberofparticle.csv");
//    myfile << Size;
//    myfile.close();
//    theTrap.extractTrapParameters("Trap Parameters.csv");
//        electronPlasma.extractPlasmaParameters("ElectronsParameters.csv");
//        theTrap.extractPlasmasHistories("");
//        electronPlasma.extractSelfPotential("SelfPotential.csv");
//        theTrap.extractTrapPotential("TrapPotential.csv");


    
//    electronPlasma.extractInitialDensitybinary("herebin.bin");
  return 0;
    

}
    //electrons.extractselfpotential11 ("./Semester 2/Potentialenergy/correctionpotential11.csv", 0.02669);
   //   electrons.extractselfpotential12 ("./Semester 2/Potentialenergy/correctionpotential12.csv", 0.02669);
  //   electrons.loadSingleRing( -100 * ePos, 0, 0.02669, 0);
  //   electrons.extracttotalpotentialenergy("./Semester 2//Potentialenergy/potentialenergycorrection.csv", 0,0.02669);
// electrons.extractPlasmaParameters("./Semester 2/Potentialenergy/ ElectronsParameters.csv");
/*
std::vector<Electrode> trapElectrodes;
   trapElectrodes.push_back(Electrode(aLength, -50));
   trapElectrodes.push_back(Electrode(aLength, 0));
   trapElectrodes.push_back(Electrode(aLength, 0));
   trapElectrodes.push_back(Electrode(aLength, -50));
   std::vector<double> trapGaps(3, aGap);
   PenningTrap theTrap(trapRadius, trapElectrodes, trapGaps, Nz, Nr);
   Plasma electronPlasma(theTrap, "electrons", massE, -ePos);
   electronPlasma.loadProfile(100,-1e-12,3.5,3.5,1e5);
for(int i = 0; i < 300; ++i)
    {
    theTrap.saveStates(i);
    theTrap.movePlasmas(18e-9);
    }

 
   theTrap.extractTrapParameters("Trap Parameters.csv");
   electronPlasma.extractPlasmaParameters("ElectronsParameters.csv");
   theTrap.extractPlasmasHistories("");
   electronPlasma.extractSelfPotential("SelfPotential.csv");
    theTrap.extractTrapPotential("TrapPotential.csv");


 for (double voltage=0; voltage< 12; ++voltage )
    yes
    /*/
 
