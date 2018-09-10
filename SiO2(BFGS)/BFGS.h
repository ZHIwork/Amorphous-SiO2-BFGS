#ifndef BFGS_H_
#define BFGS_H_

#include <iostream>
#include "CONSTANT.h"
#include "EnergyForce.h"
#include "Tool.h"

class BFGS
{
private:
	//all of these arrays are in Cartesian Coordinates
	int dimension;
	T *inverseB0;
	T *x_before, *x_after;
	T *inverseB_before, *inverseB_after;
	T *g_before, *g_after;
	T *S_before;
	T *theta, *y;
	T Alpha, Beta, Gamma, t;
	T *temp2Dimension, *tempDimension,*tempOne, *I;
public:
	BFGS();
	~BFGS();

public:
	void LineSearch();
	T GetAlpha(const T &alpha, const T*Sk, const T*gk);
	void free1Darray(T *array);
	void setZero(T *a, const int &N);
	void assignValuefor1Darray(const T *a, T *b, const int &N);
	void add(const T *a, const T *b, T *c, const int &N);
	void subtract(const T *a, const T *b, T *c, const int &N);
	void numMultiply(const T &num, const T *a, T *b, const int &N);
	T specialMultiply(const T *a, const T *b, const int &N);
	void convertDirecttoCartesianCoordinates(const T *coorDirect, T *coorCartesian, const int &N);
	void convertandLimitCarttoDirectCoordinates(const T *coorCartesian, T *coorDirect, const int &N);
	void limitCoor(T *a, const int &N);
};


#endif /* BFGS_H_ */