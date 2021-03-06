#ifndef GROUP_H_
#define GROUP_H_

#include <iostream>
#include <random>
#include "CONSTANT.h"
#include "Tool.h"

class Group
{
private:
	int orderofO;
	T directCoorO[3], cartesianCoorO[3];
	T directCoorSi[3], cartesianCoorSi[3];
	T allCoorO[27][3];
	int orderSi[2], ordernotcenterO[2][3];
public:
	Group();
	~Group();

public:
	void Graph();
	void allPossibleCoordinates(const T (&coor)[3], T (&allCoor)[27][3]);//for direct coordinate
	void groupfromGraph();
	void cutBond(int graph[][128]);
	int generator_int(int begin, int end);
	
};

#endif /* GROUP_H_ */
