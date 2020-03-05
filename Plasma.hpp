/*
Written by: Daniel Duque and Clementine Domine
Last modified on 10 Dec 2019
Declarations for the Plasma class
This file contains a corresponding source file.
*/

#ifndef PLASMA_HPP
#define PLASMA_HPP

#include"Constants.hpp"
#include"PenningTrap.hpp"
#include</Users/clementine/Documents/Documents – Clementine’s MacBook Pro/C++/eigen-eigen-323c052e1731/Eigen/SparseLU>
#include<array>
#include<vector>

class MacroRing
{
private:
    int posR;//Fixed at one of the grid divisions along r
    double posZ; //Position along z
    double speed; //speed along z, there is no r transfer.
    std::vector<double> historySpeed;
    std::vector<double> historyZ;
public:
    MacroRing(int r, double z, double aSpeed);
    ~MacroRing();
    int getR() const;
    double getZ() const;
    double getSpeed() const;
    void setZ(double newZ);
    void setSpeed(double newSpeed);
    void saveState(); //Stores the position and speed.
    void reserve(int desired); //Saves space for how many times you will save state.
    void printPositions(std::ofstream& file) const;
    void printSpeeds(std::ofstream& file) const;
};

class Plasma
{
private:
    PenningTrap& refTrap;
    const std::string name;//A name to identify the output files e.g. Electrons, Antiprotons, etc.
    const double mass; //Mass of type of particles that form the plasma. It is NOT the total mass of the plasma.
    const double charge; //Not the net charge of the plasma. It is the charge of the particle that forms the plasma.
    double chargeMacro;//Charge of each MacroRing, all rings have the same charge regardless of r position.
    double massMacro; //Mass of each MacroRing.
    std::vector<MacroRing> rings;
    //double temperature;
    Eigen::VectorXd selfPotential{ refTrap.Nz * refTrap.Nr + refTrap.Nr };
    Eigen::VectorXd RHS{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// = - charge density / permittivity of free space
    void updateRHS();
    void solvePoisson();
    void moveRings(double deltaT); //Move all particles by deltaT seconds
    void saveState(); //Store state of each of the rings
    void reserve(int desired);
    void extractHistory(std::string preName) const;
    double temperature;
    double totalCharge;
    double var;
    double var1;
    double var2;
    
    Eigen::VectorXd initialDensity1{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// Initial density with which the plasma was loaded
    Eigen::VectorXd initialDensity2{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// Initial density with which the plasma was loaded
     Eigen::VectorXd initialDensity3{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// Initial density with which the plasma was loaded
    Eigen::VectorXd initialDensity4{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// Initial density with which the plasma was loaded
    Eigen::VectorXd initialDensity5{ refTrap.Nz * refTrap.Nr + refTrap.Nr };// Initial density with which the plasma was loaded
    double weight;
public:
    Plasma(PenningTrap& trap, std::string name, double mass, double charge);
    ~Plasma();
    friend PenningTrap;
    void extractSelfPotential(std::string fileName) const;//Store the potential field in a document called filename
    void extractPlasmaParameters(std::string filename) const;//Extract the plasma parameters
    
    //Loading routines
    /*----------------------------------------------------------------------------------
    These routines below are the only thing that you should change/add to the code.
    Add a different loading routine depending on how you want the plasma to be loaded
    into the trap.
    The template for a loading routine should be:
    1) Check that input makes sense.
    2) Assign charge and mass of a MacroRing (same for all).
    3) Clear the rings vector (This is like cleaning the trap before loading a plasma).
    4) Push the MacroRings into the rings vector according to the loading you want.
    5) Solve Poisson's Equation. Just call the method; literally one line: solvePoisson();
    The code only requires the loading to populate the rings vector i.e. telling the code how
    and where are the particles initially.
    -------------------------------------------------------------------------------------*/
    void loadOneDUniform(int numMacro, double chargeMacro, double lengthLine, int r); //equally spaced particles at fixed r in a length centred in the center of the trap.
    void loadSingleRing(double chargeMacro, int r, double Z, double speed); //Single ring at a position with a velocity
    std::vector<double> getselfpotential11(double Z); //get the Selfpotencial of one particle at all the R ( not Z needs to be on a grid point)
       std::vector<double> getselfpotential12(double Z); //get the Selfpotencial of one particle at all the R ( not Z needs to be on a grid point)
    void extractselfpotential11 (std::string fileName, double Z);
    void extractselfpotential12 (std::string fileName, double Z);
    std::vector<double> gettotalpotentialenergy(int r,double Z);// all the double Z in this cases should be the same
    void extracttotalpotentialenergy(std::string fileName, int r,double Z);
    void loadProfile(double aTotalcharge, double b, double n,double T);
    void extractInitialDensity(std::string fileName) const;
    void estimateDensityProportions();
    void fitDensityProportionToProfile( double n, double b);
    void normalizeDensityToTotalCharge();
    void estimateDensityProportions2();
    void fitDensityProportionToProfile2( double n, double b);
    void normalizeDensityToTotalCharge2();
    void extractInitialDensitybinary(std::string fileName) const;
    
    };
#endif
