//============================================================================
// Name        : main.cpp
// Author      : Zhi Song
// Version     : 1.0.0
// Copyright   : Your copyright notice
// Description : Create generic code to generate amorphous system base on some
//				 crystal structure like SiO2 in C++
//============================================================================

#include <iostream>
#include <chrono>
#include "CONSTANT.h"
#include "IO.h"
#include "Group.h"
#include "EnergyForce.h"
#include "BFGS.h"

using namespace std;

void freeAllMemory();


int main() {

	//cout << "start" << endl;

	IO mainIO;
	mainIO.readfromFile();
	
	Group mainGroup;
	mainGroup.Graph();

	EnergyForce mainEF;
	Tool mainTool;
	
	//cut bonds using graph
	//mainGroup.cutBond(graph);
	//group from graph
	mainGroup.groupfromGraph();
	//mainEF.obtainEnergyForce(atomCoordinate, atomAcceleration, totalEnergy);
	//cout <<" Energy: "<<totalEnergy<<endl;
	

	for (int i = 0; i < cutBondsTimes; ++i) {

		//auto start = chrono::steady_clock::now();

		//store parameters
		mainTool.copyGraph(graph, oldGraph, atomQuantity[1]);//store graph to old graph
		mainTool.copyCoordinates(atomCoordinate, atomOldCoordinate, totalAtomQuantity * 3);//store atom coordinates to old atom coordiantes
		oldEnergy = totalEnergy;//store total energy to old energy


		//cut bonds using graph
		mainGroup.cutBond(graph);

		//group from graph
		mainGroup.groupfromGraph();

		//do linesearch (BFGS)
		BFGS mainBFGS;
		mainBFGS.LineSearch();

		//MC determination 
		auto p = mainTool.MC_probability(oldEnergy, totalEnergy);	

		cout << "Old energy: " << oldEnergy << "    New Energy: " << totalEnergy << endl;
	
		//if p==0, reuse old energy, old graph and old coordinates
		mainTool.implementMC(p);

		//output file
		mainIO.output(i);

	}

	//free memory
	freeAllMemory();


	return 0;
}



void freeAllMemory()
{
	//free memory of latticeVector
	freeLatticeVector();

	//free memory of inverseLatticeVector
	freeInverseLatticeVector();

	//free memory of atomCoordinate
	freeAtomCoordinate();

	//free memory of atomOldCoordinate
	freeAtomOldCoordinate();

	//free memory of atomAcceleration
	freeAtomAcceleration();
}

