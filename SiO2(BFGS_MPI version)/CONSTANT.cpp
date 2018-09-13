#include "CONSTANT.h"

//IO interface
const std::string importFileName = "SiO2.vasp";//input file name.
const std::string exportFileName = "final.vasp";//output file name.


//simulation parameters
const int converganceTimes = 100;//when minimize energy, the number of times to converge.
const int cutBondsTimes = 1000;//when simulate amorphous system, the number of times to cut bonds.


//Environment parameters
const T temperature = 6000.0;//temperature that used to form amorphous system.


//Constant that uses in equations
const T pi = 3.141592653589793;//pi constant

const T BoltzmannConstant = 8.617330350e-5;//Boltzmann Constant

const T massSi = 46.82e-27;//mass of Si
const T massO = 26.67e-27;//mass of O

const T b0 = 1.55148;
const T kb = 26.96;
const T kthetaSi = 1.685;
const T kthetaO = 0.58;
const T thetaSi = 109.47;
const T thetaO = 180.0;
const T cos0_Si = cos(thetaSi * pi / 180.0);
const T cos0_O = cos(thetaO * pi / 180.0);
const T repulsiveCoeforSiSi = 8.0;
const T repulsiveCoeforOO = 2.0;


//POSCAR CONSTANT
T latticeConstant = 0.0; //universal scaling factor (lattice constant),
						 //which is used to scale all lattice vectors and all atomic coordinates. 
T latticeVector[3 * 3]{ 0 }; //three lattice vectors defining the unit cell of the system are given

T inverseLatticeVector[3 * 3]{ 0 };//inverse lattice vector

const int atomType = 2; //the quantity of atom's types (need to be correct when use different POSCAR)
std::string atomName[atomType]; //all atom's names from POSCARS
const int atomQuantity[atomType]={128, 64}; //the quantity of each type of atoms (O, Si orders)
const int totalAtomQuantity = 192; // the total number atoms (need to be correct when use different POSCAR)

T atomCoordinate[totalAtomQuantity * 3]{ 0 };//store all coordinates of x, y, z into this dynamic array. 

T atomOldCoordinate[totalAtomQuantity * 3]{ 0 };//store all coordinates of x, y, z into this dynamic array. 

T atomAcceleration[totalAtomQuantity * 3]{ 0 };//store all acceleration coordinates of x, y, z into this dynamic array. 

int startPoint[atomType] = {0, 128*3};//start point for different type atoms, need to be changed when use different POSCAR


//For Group storing
int graph[64][128]{ 0 };//first one is the number of Si and P, second is the number of O.
int oldGraph[64][128]{ 0 };//first one is the number of Si and P, second is the number of O.
int groupSi[64][4];//array to store all bonded O with each Si 
int groupO[128][2];//array to store all bonded Si with each O


//
T totalEnergy = 0;//total energy of whole system
const T frequency = 1e-15;
T oldEnergy = 0;
T infiniteNumber = 1e8;
