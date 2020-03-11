/*
Written by: Daniel Duque and Clementine Domine
Last modified on 10 Dec 2019
Definitions for the Plasma class
This file contains a corresponding source file.
*/
#include"Plasma.hpp"
#include<iostream>
#include<fstream>

MacroRing::MacroRing(int r, double z, double aSpeed) : posR(r), posZ(z), speed(aSpeed)
{
}
int MacroRing::getR() const
{
    return posR;
}
double MacroRing::getZ() const
{
    return posZ;
}
double MacroRing::getSpeed() const
{
    return speed;
}
void MacroRing::setZ(double newZ)
{
    posZ = newZ;
}
void MacroRing::setSpeed(double newSpeed)
{
    speed = newSpeed;
}
void MacroRing::saveState()
{
    historyZ.push_back(posZ);
    historySpeed.push_back(speed);
}
void MacroRing::reserve(int desired)
{
    historyZ.reserve(desired);
    historySpeed.reserve(desired);
}
void MacroRing::printPositions(std::ofstream& file) const
{
    if (!historyZ.empty())
    {
        file << posR;
        file << ",";
        for (unsigned int i = 0; i < historyZ.size() - 1; ++i)
        {
            file << historyZ[i] << ",";
        }
        file << historyZ.back() << "\n";
    }
}
void MacroRing::printSpeeds(std::ofstream& file) const
{
    for (unsigned int i = 0; i < historySpeed.size(); ++i)
    {
        file << historySpeed[i];
        i < historySpeed.size() - 1 ? file << "," : file << "\n";
    }
}

/*----------------------------------------------------------------------------------------------------------
Plasma class
----------------------------------------------------------------------------------------------------------*/

Plasma::Plasma(PenningTrap& trap, std::string aName, double aMass, double aCharge)
    : refTrap(trap) , name(aName), mass(aMass), charge(aCharge)
{
    refTrap.addPlasma(*this);
}
Plasma::~Plasma()
{
    //Should remove itself from the plasma, look later into this
}
void Plasma::updateRHS()
{
    //The RHS is: (minus) chargeDensity/epsilon
    //Performs a First-order weighting (Area weighting)
    RHS.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    double hr{ refTrap.hr };
    double hz{ refTrap.hz };
    for (const MacroRing& aRing : rings)
    {
        int indexR{ aRing.getR() };
        int indexZ{ (int)floor(aRing.getZ() / hz) };
        int indexRHS{ (refTrap.Nz + 1) * indexR + indexZ };
        double z{ aRing.getZ() - indexZ * hz };//Distance from left grid point to particle
        double weightFactor{ z / hz };//Weight going to grid point on the right (number between 0 and 1)
        RHS.coeffRef(indexRHS) += -macroChargeDensity * (1 - weightFactor) / epsilon;//Point to the left
        RHS.coeffRef(indexRHS + 1) += -macroChargeDensity * weightFactor / epsilon;//Point to the right
    }
}
void Plasma::solvePoisson()
{
    updateRHS();
    selfPotential = refTrap.solver.solve(RHS);
}
void Plasma::moveRings(double deltaT)
{
    for (unsigned int i = 0; i < rings.size(); )
    {
        //Leapfrog method
        
        
        
        double vNew{ deltaT * refTrap.getEField(rings[i].getR(), rings[i].getZ()) * charge / mass + rings[i].getSpeed() };
        double zNew{ deltaT * vNew + rings[i].getZ() };
        //remove elements than escape the trap
        if (zNew < refTrap.getLength() && zNew > 0)
        {
            rings[i].setZ(zNew);
            rings[i].setSpeed(vNew);
            ++i;
        }
        else
        {
            std::swap(rings[i], rings.back());//check that this std::swap works efficiently as expected
            rings.pop_back();
        }
    }
}
void Plasma::extractSelfPotential(std::string fileName) const
{
    Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
    std::ofstream newFile;
    newFile.open(fileName);
    newFile << selfPotential.format(fastFullPrecision);
    newFile.close();
}
void Plasma::extractPlasmaParameters(std::string fileName) const
{
    std::ofstream newFile;
    newFile.open(fileName);
    newFile << mass << '\n';
    newFile << charge << '\n';
    newFile << chargeMacro;
    newFile.close();
}
void Plasma::extractHistory(std::string preName) const
{
    std::ofstream newPositions, newSpeeds;
    newPositions.open(preName + "Positions" + name + ".csv");
    newSpeeds.open(preName + "Speeds" + name + ".csv");
    for (const MacroRing& aRing : rings)
    {
        aRing.printPositions(newPositions);
        aRing.printSpeeds(newSpeeds);
    }
    newPositions.close();
    newSpeeds.close();
}

double Plasma::getPotentialEnergy() const
{
    double potentialEnergy{ 0 };
    for (const MacroRing& aRing : rings)
    {
        potentialEnergy += refTrap.getTotalPhi(aRing.getR(), aRing.getZ()) * (aRing.getR() == 0 ? chargeMacro : aRing.getR() * 8 * chargeMacro);
    }
    return potentialEnergy / 2;
}
void Plasma::saveState()
{
    for (MacroRing& aRing : rings)
    {
        aRing.saveState();
    }
}
void Plasma::reserve(int desired)
{
    for (MacroRing& aRing : rings)
    {
        aRing.reserve(desired);
    }
}
void Plasma::loadOneDUniform(int numMacro, double aTotalCharge, double lengthLine, int r)
{
    //Creates equally spaced charges at r = 0 along a centred line of length lengthLine
    if (lengthLine >= refTrap.getLength() || r >= refTrap.Nr - 1)
    {
        throw std::logic_error("Length of charge must be less than the length of the trap and r inside the trap");
    }
    if (aTotalCharge * charge < 0)
    {
        throw std::logic_error("Charge of MacroRing and the plasma type must have the same sign");
    }
    aTotalCharge /= numMacro; //Charge of a single macro-ring
    chargeMacro = r == 0 ? aTotalCharge : aTotalCharge / (8 * r);
    macroChargeDensity = 4 * chargeMacro / (PI * refTrap.hz * refTrap.hr * refTrap.hr);
    rings.clear();
    rings.reserve(numMacro);
    //Divide lengthLine in numMacro + 1 cells, which corresponds to numMacro + 2 points
    //Put the particles along the points except for the first and last point
    double start{ (refTrap.getLength() - lengthLine) / 2 };
    double hz{ lengthLine / (numMacro + 1) };
    for (int i = 1; i <= numMacro; ++i)
    {
        rings.push_back(MacroRing(r, start + i * hz, 0));
    }
    solvePoisson();
}
void Plasma::loadSingleRing(double aTotalCharge, int r, double Z, double speed)
{
    if (Z >= refTrap.getLength() || Z <= 0 || r >= refTrap.Nr - 1)
    {
        throw std::logic_error("Input r,z is not inside the trap");
    }
    if (aTotalCharge * charge < 0)
    {
        throw std::logic_error("Charge of MacroRing and the plasma type must have the same sign");
    }
    chargeMacro = r == 0 ? aTotalCharge : aTotalCharge / (8 * r);
    macroChargeDensity = 4 * chargeMacro / (PI * refTrap.hz * refTrap.hr * refTrap.hr);
    rings.clear();
    rings.push_back(MacroRing(r, Z, speed));
    solvePoisson();
}

std::vector<double> Plasma::getselfpotential11(double Z)
{   std::vector<double> SPE;
    for (int i = 0; i < refTrap.Nr-1; ++i)
    {   int indexZ{ (int)floor(Z/ refTrap.hz) };
        loadSingleRing( -100 * ePos, i, indexZ*refTrap.hz, 0);
        double spe= selfPotential(indexZ+(i*(refTrap.Nz+1)));
        SPE.push_back(spe);
    }
    return SPE;
}

std::vector<double> Plasma::getselfpotential12(double Z)
{   std::vector<double> SPE12;
    for (int i = 0; i < refTrap.Nr-1; ++i)
    {int indexZ{ (int)floor(Z/ refTrap.hz) };
        loadSingleRing( -100 * ePos, i, indexZ*refTrap.hz, 0);
        double spe= selfPotential((indexZ+1)+(i*(refTrap.Nz+1)));
        SPE12.push_back(spe);
    }
    return SPE12;
}


void Plasma::extractselfpotential11 (std::string fileName, double Z)
{
    std::ofstream newFile;
    newFile.open(fileName);
    std::vector<double>  a;
    a = getselfpotential11(Z);
    
    for (unsigned int i = 0; i < a.size(); ++i)
        {
            newFile << a[i];
            i < a.size() - 1 ? newFile << "," : newFile << "\n";
        }
    newFile.close();
}

void Plasma::extractselfpotential12 (std::string fileName, double Z)
{
    std::ofstream newFile;
    newFile.open(fileName);
    std::vector<double>  a;
    a = getselfpotential12(Z);
    
    for (unsigned int i = 0; i < a.size(); ++i)
        {
            newFile << a[i];
            i < a.size() - 1 ? newFile << "," : newFile << "\n";
        }
    newFile.close();
}

std::vector<double> Plasma::gettotalpotentialenergy(int r,double Z){
std::vector<double> Potentialenergy;
        double ps11 = getselfpotential11(Z)[r];
        double ps12 = getselfpotential12(Z)[r];
        loadOneDUniform(10000, -100 * ePos, 0.01, 0);
    for (const MacroRing& aRing : rings)
        {
        int indexZ{ (int)floor(aRing.getZ() / refTrap.hz) };
        double dz{ aRing.getZ() - indexZ * refTrap.hz };
        double weightFactor{ dz / refTrap.hz };
        double ps = (pow((1 - weightFactor),2)+pow((weightFactor),2))*ps11+2*((1 - weightFactor)* weightFactor * ps12);
            
        double  phiselfleft = selfPotential(indexZ+(r*(refTrap.Nz+1)));
        double  phiselftright = selfPotential(indexZ+1+(r*(refTrap.Nz+1)));
        double  phiT= phiselfleft*(1 - weightFactor)+phiselftright*(weightFactor);//Point to the left
        
            
        double phitrapleft =refTrap.potentialsVector(indexZ+(r*(refTrap.Nz+1)));
        double phitrapright =refTrap.potentialsVector(indexZ+1+(r*(refTrap.Nz+1)));
        double phiTraptot= phitrapleft*(1 - weightFactor)+phitrapright*(weightFactor);//Point to the left
        
        double
        specificCharge{r== 0 ? chargeMacro : 8 * chargeMacro * r};
        double phireal= specificCharge*(ps);
        Potentialenergy.push_back(phireal);
                                                        
        }
    return Potentialenergy;
}

void Plasma::extracttotalpotentialenergy(std::string fileName, int r,double Z)
 {
    std::ofstream newFile;
    newFile.open(fileName);
    std::vector<double>  b;
    b = gettotalpotentialenergy(r,Z);
         for (unsigned int i = 0; i < b.size(); ++i)
             {
                 newFile << b[i];
                 i < b.size() - 1 ? newFile << "," : newFile << "\n";
             }
         newFile.close();
     }


void Plasma::estimateDensityProportions()
{
    int pointsZ = refTrap.Nz + 1;
    int pointsR = refTrap.Nr;
     //For each grid point, estimate density from Total Potential
     for (int indexR = 0; indexR < pointsR; ++indexR)
     {
         double phiCentralR{ refTrap.getTotalPhi(indexR, refTrap.lengthTrap / 2) };
         for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
         {
             initialDensity2.coeffRef(pointsZ * indexR + indexZ) = exp(-(charge / (KB * temperature)) * (refTrap.getTotalPhi(indexR, indexZ) - phiCentralR));
         }
     }
}

void Plasma::fitDensityProportionToProfile( double n, double b)
{
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
         {
             double initialDensity2normvalue {0};
             
             for (int indexZ = 0; indexZ < (refTrap.Nz + 1); ++indexZ)
             {
             //Convert from volume density into flat space density
                 initialDensity2normvalue += initialDensity2.coeffRef((refTrap.Nz +1)* indexR + indexZ) * refTrap.hz;
             }
             //this number becomes very big
             double Factornorm2 = ((exp(-pow((indexR* refTrap.hr *1000 /b),n))/initialDensity2normvalue));
             //This 1000 is because the profile is given in milimeters
             
             for (int indexZ = 0; indexZ < (refTrap.Nz + 1 ); ++indexZ)
             {
                 initialDensity2.coeffRef((refTrap.Nz + 1 )* indexR + indexZ) *= Factornorm2;
             }
                  
         }
    
}

void Plasma::normalizeDensityToTotalCharge()
{   double volume ;
    double currentCharge{ 0 };
       for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
       {
           //First calculate what the total charge is right now
           if (indexR== 0)//Volume of the MacroRing
           {
              volume = PI * refTrap.hz * refTrap.hr * refTrap.hr / 4;
           }
           else
           {
                volume =refTrap.hz * refTrap.hr* 2 * PI * indexR *  refTrap.hr;
           }
           for (int indexZ = 0; indexZ < (refTrap.Nz + 1 ); ++indexZ)
           {
               currentCharge += (initialDensity2.coeffRef((refTrap.Nz + 1)* indexR + indexZ) * volume);
           }
        
       }
       //Now multiply every density by the appropriate correction to get the expected total charge
       double correction= totalCharge/currentCharge;
       for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
       {
           for (int indexZ = 0; indexZ < refTrap.Nz+1; ++indexZ)
           {
               initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ) *= correction;
           }
    }
    
}

void Plasma::estimateDensityProportions2()
{
    int pointsZ = refTrap.Nz + 1;
    int pointsR = refTrap.Nr;
     //For each grid point, estimate density from Total Potential
     for (int indexR = 0; indexR < pointsR; ++indexR)
     {
         double phiCentralR{ refTrap.getTotalPhi(indexR, refTrap.lengthTrap / 2) };
         for (int indexZ = 0; indexZ < pointsZ; ++indexZ)
         {
             initialDensity4.coeffRef(pointsZ * indexR + indexZ) = exp(-(charge / (KB * temperature)) * (refTrap.getTotalPhi(indexR, indexZ) - phiCentralR));
         }
     }
}

void Plasma::fitDensityProportionToProfile2( double n, double b)
{
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
         {
             double initialDensity4normvalue {0};
             
             for (int indexZ = 0; indexZ < (refTrap.Nz + 1); ++indexZ)
             {
             //Convert from volume density into flat space density
                 initialDensity4normvalue += initialDensity4.coeffRef((refTrap.Nz +1)* indexR + indexZ) * refTrap.hz;
             }
             //this number becomes very big
             double Factornorm2 = ((exp(-pow((indexR* refTrap.hr *1000 /b),n))/initialDensity4normvalue));
             //This 1000 is because the profile is given in milimeters
             
             for (int indexZ = 0; indexZ < (refTrap.Nz + 1 ); ++indexZ)
             {
                 initialDensity4.coeffRef((refTrap.Nz + 1 )* indexR + indexZ) *= Factornorm2;
             }
                  
         }
    
}

void Plasma::normalizeDensityToTotalCharge2()
{   double volume ;
    double currentCharge{ 0 };
       for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
       {
           //First calculate what the total charge is right now
           if (indexR== 0)//Volume of the MacroRing
           {
              volume = PI * refTrap.hz * refTrap.hr * refTrap.hr / 4;
           }
           else
           {
                volume =refTrap.hz * refTrap.hr* 2 * PI * indexR *  refTrap.hr;
           }
           for (int indexZ = 0; indexZ < (refTrap.Nz + 1 ); ++indexZ)
           {
               currentCharge += (initialDensity4.coeffRef((refTrap.Nz + 1)* indexR + indexZ) * volume);
            
           }
        
       }
       //Now multiply every density by the appropriate correction to get the expected total charge
       double correction= totalCharge/currentCharge;
       for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
       {
           for (int indexZ = 0; indexZ < refTrap.Nz+1; ++indexZ)
           {
               initialDensity4.coeffRef((refTrap.Nz+1)* indexR + indexZ) *= correction;
           }
    }
    
}

//the variance is not used to decide th enumber of iteration in the loop
//n0 density at the center, N number of particles , a b and n paramter of the radial fit.
void Plasma::loadProfile(double aTemperature, double aTotalCharge, double n, double b, int numMacro)
{
    rings.clear();
    //Check input makes sense
    if (aTemperature <= 0|| n <= 0 || b <= 0)
    {
        throw std::logic_error("Temperature, shape, and scale all need to be positive");
    }
    if (aTotalCharge *charge < 0)
    {
        throw std::logic_error("Total charge and the plasma type charge must have the same sign");
    }
    temperature = aTemperature;
    totalCharge = aTotalCharge;
    
    //Declaration and Initialisation of the variables
    double C=0.997;
    
     // Factor for normalisation for th radial profil with n0 in center
    initialDensity1.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    initialDensity2.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    initialDensity3.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    initialDensity4.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    initialDensity5.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);

    //Solve Poisson's equation for the charge density
    //Then estimate a thermal equilibrium density from that potential (normalized to totalCharge). Should add a double return here, to get how much it changed
    //Tune the obtained density according to the expected profile
    //Repeat until it is both self consistent and in agreement with expected profile
    //Change this for a do while loop after you have figured out what conditions to impose to decide wether to repeat or not
    // Start iteration procedure
var=1e6;
var1=0;
var2=0;
double x2;
double x;
//int numIterations{10000};
//for (int i = 0; i < numIterations; ++i)
do{
    //Format the density as the right hand side of Poissons equation
    for (int j = 0; j < refTrap.Nz * refTrap.Nr + refTrap.Nr; ++j)
        {
          initialDensity1.coeffRef(j) = -initialDensity1.coeffRef(j) / epsilon;
        }
        //Solve Poisson's equation using Denstiy1
    selfPotential = refTrap.solver.solve(initialDensity1);
    //first estimate the proportions based on thermal equilibrium assumption
    estimateDensityProportions();
    //Then weight them appropriately to match the given profile
    fitDensityProportionToProfile(n, b);
    normalizeDensityToTotalCharge();
   
    //Calculate the new estimate as the linear combination of the orginal density1 and the first estimate of the density2
    x=0;
    x2=0;
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
    {
        for (int indexZ = 0; indexZ < refTrap.Nz+1; ++indexZ)
        {
            x2 += pow((-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ))-initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ),2)*indexR;
            x += pow((-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ)),2)+(pow(initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ),2))*indexR;
        
            initialDensity3.coeffRef((refTrap.Nz+1)* indexR + indexZ) =((C)*(-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ)))+(1-C)*initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ);
        }
      }
    
    var1= x2/x;
    
    for (int j = 0; j < refTrap.Nz * refTrap.Nr + refTrap.Nr; ++j)
        {
          initialDensity3.coeffRef(j) = -initialDensity3.coeffRef(j) / epsilon;
        }
    //Solve Poisson's equation using Denstiy1
    selfPotential = refTrap.solver.solve(initialDensity2);
    estimateDensityProportions2();
    //Then weight them appropriately to match the given profile
    fitDensityProportionToProfile2(n, b);
    normalizeDensityToTotalCharge2();
    
    x=0;
    x2=0;
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
    {
        for (int indexZ = 0; indexZ < refTrap.Nz+1; ++indexZ)
        {
            x2 += pow((-epsilon*initialDensity3.coeffRef((refTrap.Nz+1)* indexR + indexZ))-initialDensity4.coeffRef((refTrap.Nz+1)* indexR + indexZ),2)*indexR;
            x += pow((-epsilon*initialDensity3.coeffRef((refTrap.Nz+1)* indexR + indexZ)),2)+pow(initialDensity4.coeffRef((refTrap.Nz+1)* indexR + indexZ),2)*indexR;

        }
    }
    var2= x2/x;
       
    if (var2>var1)
    {
        C=C*0.99999999997;
    }
    if(var2<var1)
    {
        C=C*1.0000000003;
    }
        
    x=0;
    x2=0;
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
    {
       for (int indexZ = 0; indexZ < refTrap.Nz+1; ++indexZ)
       {
           x2 += pow((-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ))-initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ),2)*indexR;
           x += pow((-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ)),2)+(pow(initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ),2))*indexR;
       
           initialDensity5.coeffRef((refTrap.Nz+1)* indexR + indexZ) =((C)*(-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ)))+(1-C)*initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ);
       }
     }

    extractInitialDensity("Initial Density Estimate.csv");
    var=x2/x;
    for (int j = 0; j < refTrap.Nz * refTrap.Nr + refTrap.Nr; ++j)
          {
            initialDensity1.coeffRef(j) = initialDensity5.coeffRef(j);
          }
 } while (var>1e-2);
//Ok, now I have a density grid. I need now to populate macro-particles to match this grid density.
    //How many particles am I placing at each radial index

std::vector<std::vector<double>> cumulativeAtR;
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
    {
        std::vector<double> cumulative;
        double volume;
        if (indexR == 0)
        {
            volume = PI * refTrap.hz *  refTrap.hr *  refTrap.hr / 4;
        }
        else
        {
            volume = refTrap.hz *  refTrap.hr * 2 * PI * indexR *  refTrap.hr;
        }
        double aChargeR{ 0 };
        for (int indexZ = 0; indexZ < refTrap.Nz + 1; ++indexZ)
        {
            aChargeR += volume * initialDensity5.coeffRef((refTrap.Nz + 1) * indexR + indexZ);
            cumulative.push_back(aChargeR);
        }
        cumulativeAtR.push_back(cumulative);
    }
    double futureMacroCharge{ cumulativeAtR[0].back() };
    for (unsigned int i = 1; i < cumulativeAtR.size(); ++i)
    {
        futureMacroCharge += cumulativeAtR[i].back() / (8 * i);
    }
    chargeMacro = futureMacroCharge / numMacro;
    macroChargeDensity = 4 * chargeMacro / (PI * refTrap.hz * refTrap.hr * refTrap.hr);
    std::vector<int> numAtR;
    numAtR.push_back((int)round(cumulativeAtR[0].back() / chargeMacro));
    for (unsigned int i = 1; i < cumulativeAtR.size(); ++i)
    {
        numAtR.push_back((int)round(cumulativeAtR[i].back() / (8 * i * chargeMacro)));
    }
    //Now I know how many particles I want to place at each radius.
    //Distribute them using an inverse transform sampling
    //But without any randomness, just equally space the uniform distribution
    rings.clear();
    rings.reserve(numMacro);//this is not exactly the number of particles we will load, but it is close
    //Use std::random to produce boltzmann distributed speeds
    //this are simply gaussian with the appropriate std deviation as we are looking at 1 dimension only
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, sqrt(KB * aTemperature / mass));
    for (int indexR = 0; indexR < refTrap.Nr; ++indexR)
    {
        double deltaQ{ cumulativeAtR[indexR].back() / (numAtR[indexR] + 1) };
        int currentIndex{ 0 };
        for (int i = 0; i < numAtR[indexR]; ++i)
        {
            double currentInvert{ deltaQ * (i + 1) };
            while (abs(cumulativeAtR[indexR][currentIndex]) < abs(currentInvert))
            {
                currentIndex++;
            }
            //The point to invert is between currentIndex and currentIndex - 1
            //Assume linear interpolation between this two points
            double invertedPos{ (currentIndex - 1) * refTrap.hz + refTrap.hz / 2 + refTrap.hz * (currentInvert - cumulativeAtR[indexR][currentIndex - 1]) / (cumulativeAtR[indexR][currentIndex] - cumulativeAtR[indexR][currentIndex - 1]) };
            rings.push_back(MacroRing(indexR, invertedPos, distribution(generator)));
        }
    }
    std::cout << "Loading " << rings.size() << " macro-particles from which " << numAtR[0] << " are at r=0.\n";
    
    solvePoisson();
}

void Plasma::extractInitialDensity(std::string fileName) const
    {
            Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
            std::ofstream newFile;
            newFile.open(fileName);
            newFile << initialDensity5.format(fastFullPrecision);
            newFile.close();
                         }


double Plasma::extractsizeringvector()
{   double size = rings.size();
    return size;
}

void Plasma::extractInitialDensitybinary(std::string fileName)
{
    loadOneDUniform(10000, -100 * ePos, 0.01, 0);
//--------------------------------------------------------------------------------------------------------------------------------------------------
    std::ofstream newFile;
    newFile.open(fileName, std::ios::out |  std::ios::binary);
    
    if(!newFile)
    {
        throw std::logic_error( "Cannot open file!");
    }
    for (const MacroRing& aRing : rings)
     {
         double   pos1= aRing.getZ();
         newFile.write(reinterpret_cast<char*>(&pos1),sizeof(pos1));
     }
//  double pos = 100;
//  char  C='C';
//  newFile.write((char*)(&C),sizeof(C));
//  newFile.write(reinterpret_cast<char*>(&pos),sizeof(double));
    newFile.close();
    if(!newFile.good())
    {
        throw std::logic_error("Error occurred at writing time!");
    }
//--------------------------------------------------------------------------------------------------------------------------------------------------
   std:: ifstream Readfile(fileName, std::ios::out | std::ios::binary);
    if(!Readfile)
    {
      throw std::logic_error( "Cannot open file!");
      
    }
    
    
  for (const MacroRing& aRing : rings)   {
         double   posr;
         
        Readfile.read(reinterpret_cast<char*>(&posr),sizeof(double));
        std::cout << posr;
   }
//    char Cr;
//    Readfile.read((char*)(&Cr),sizeof(Cr));
//    std::cout << Cr;
    Readfile.close();
    
    if(!Readfile.good())
    {
        throw std::logic_error("Error occurred at Readinf time!");
   }
}

//    double k = 100;
//    char buffer[100];
//    newFile.write (buffer, 100);
//    newFile.write(reinterpret_cast<char*>(&k),sizeof(k));
//z
 
