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
MacroRing::~MacroRing()
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
        double volume; //Volume of the MacroRing
        if (indexR == 0)
        {
            volume = PI * hz * hr * hr / 4;
        }
        else
        {
            volume = hz * hr * 2 * PI * indexR * hr;
        }
        double specificCharge{aRing.getR() == 0 ? chargeMacro : 8 * chargeMacro * aRing.getR()};
        
        RHS.coeffRef(indexRHS) += -specificCharge * (1 - weightFactor) / (volume * epsilon);//Point to the left
        RHS.coeffRef(indexRHS + 1) += -specificCharge* weightFactor / (volume * epsilon);//Point to the right
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
        
        double specificCharge{rings[i].getR() == 0 ? chargeMacro : 8 * chargeMacro * rings[i].getR()};
        double specificMass{rings[i].getR() == 0 ? massMacro : 8 * massMacro * rings[i].getR()};
        
        double forceOld{ specificCharge * refTrap.getEField(rings[i].getR(), rings[i].getZ()) };
        double vNew{ deltaT * forceOld / specificMass + rings[i].getSpeed() };
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
    newFile << massMacro << '\n';
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
void Plasma::loadOneDUniform(int numMacro, double aChargeMacro, double lengthLine, int r)
{
    //Creates equally spaced charges at r = 0 along a centred line of length lengthLine
    if (lengthLine >= refTrap.getLength() || r >= refTrap.Nr - 1)
    {
        throw std::logic_error("Length of charge must be less than the length of the trap and r inside the trap");
    }
    if (aChargeMacro * charge < 0)
    {
        throw std::logic_error("Charge of MacroRing and the plasma type must have the same sign");
    }
    
    chargeMacro = aChargeMacro;
    massMacro = mass *  chargeMacro / charge;
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
void Plasma::loadSingleRing(double aChargeMacro, int r, double Z, double speed)
{
    if (Z >= refTrap.getLength() || Z <= 0 || r >= refTrap.Nr - 1)
    {
        throw std::logic_error("Input r,z is not inside the trap");
    }
    if (aChargeMacro * charge < 0)
    {
        throw std::logic_error("Charge of MacroRing and the plasma type must have the same sign");
    }

    massMacro = mass * chargeMacro/ charge;
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
        
        double specificCharge{r== 0 ? chargeMacro : 8 * chargeMacro * r};
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
void Plasma::loadProfile( double aTotalcharge,  double b, double n,double aTemperature)
{
    //Check input makes sense
    if (aTemperature <= 0|| n <= 0 || b <= 0)
    {
        throw std::logic_error("Temperature, shape, and scale all need to be positive");
    }
    if (aTotalcharge *charge < 0)
    {
        throw std::logic_error("Total charge and the plasma type charge must have the same sign");
    }
    temperature = aTemperature;
    totalCharge = aTotalcharge;
    
    //Declaration and Initialisation of the variables
    double C=0.997;
    
     // Factor for normalisation for th radial profil with n0 in center
    initialDensity1.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    initialDensity2.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);
    initialDensity3.setZero(refTrap.Nz * refTrap.Nr + refTrap.Nr);

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
    selfPotential = refTrap.solver.solve(initialDensity3);
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
//               initialDensity3.coeffRef((refTrap.Nz+1)* indexR + indexZ) =((C)*(-epsilon*initialDensity1.coeffRef((refTrap.Nz+1)* indexR + indexZ)))
//                    +(1-C)*initialDensity2.coeffRef((refTrap.Nz+1)* indexR + indexZ);
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

    extractInitialDensitybinary("Daniel.dat");
    var=x2/x;
    for (int j = 0; j < refTrap.Nz * refTrap.Nr + refTrap.Nr; ++j)
          {
            initialDensity1.coeffRef(j) = initialDensity5.coeffRef(j);
          }
 } while (var>1e-1);
}

void Plasma::extractInitialDensity(std::string fileName) const
    {
            Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
            std::ofstream newFile;
            newFile.open(fileName);
            newFile << initialDensity5.format(fastFullPrecision);
            newFile.close();
                         }
void Plasma::extractInitialDensitybinary(std::string fileName) const
{
    Eigen::IOFormat fastFullPrecision(Eigen::FullPrecision, Eigen::DontAlignCols, "", "\n", "", "", "", "");
    std::ofstream newFile(fileName, std::ios::out | std::ios::binary);
    if(!newFile)
    {
        throw std::logic_error( "Cannot open file!");
    }
    newFile.open(fileName);
    for(int i = 0; i < 3; i++)
    for (int j = 0; j < refTrap.Nz * refTrap.Nr + refTrap.Nr; ++j)
                 {
    newFile.write((char *) &initialDensity5.coeffRef(j), sizeof(initialDensity5.coeffRef(j)));  //not sure what I want to write there
                  }
    newFile.close();
    
    if(!newFile.good())
    {
    throw std::logic_error("Error occurred at writing time!");
    }
    
}

