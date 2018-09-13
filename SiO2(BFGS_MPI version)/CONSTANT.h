#ifndef CONSTANT_H_
#define CONSTANT_H_
#include <string>
#include <map>
#include <cmath>

//define type
typedef double T;

//IO interface
extern const std::string importFileName;
extern const std::string exportFileName;


//simulation parameters
extern const int converganceTimes;
extern const int cutBondsTimes;


//Environment parameters
extern const T temperature;


//Constant that uses in equations
extern const T pi;

extern const T BoltzmannConstant;

extern const T massSi;
extern const T massO;

extern const T b0;
extern const T kb;
extern const T kthetaSi;
extern const T kthetaO;
extern const T thetaSi;
extern const T thetaO;
extern const T cos0_Si;
extern const T cos0_O;
extern const T repulsiveCoeforSiSi;
extern const T repulsiveCoeforOO;


//POSCAR CONSTANT
extern T latticeVector[];
extern T inverseLatticeVector[];
extern T latticeConstant;
extern const int atomType;
extern std::string atomName[];
extern const int atomQuantity[];
extern const int totalAtomQuantity;
extern T atomCoordinate[];
extern T atomOldCoordinate[];
extern T atomAcceleration[];
extern int startPoint[];


//For Group storing
extern int graph[][128];//number 128 need to be changed when changed POSCAR(also in Group.cpp need to be revised)
extern int oldGraph[][128];
extern int groupSi[][4];
extern int groupO[][2];


//
extern T totalEnergy;
extern const T frequency;
extern T oldEnergy;
extern T infiniteNumber;

#endif /* CONSTANT_H_ */





