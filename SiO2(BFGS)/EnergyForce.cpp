#include "EnergyForce.h"


EnergyForce::EnergyForce()
{
	//std::cout << "EnergyForce constructor" << std::endl;
}


EnergyForce::~EnergyForce()
{
	//std::cout << "EnergyForce deconstructor" << std::endl;
}


void EnergyForce::obtainEnergyForce(const T *coordinates, T *acceleration, T& total_energy)
{		
	//call Tool.h
	Tool EFTool;

	//initialize energy parameters
	angular_energy_Si = 0;
	angular_energy_O = 0;
	radial_energy = 0;
	repulsive_energy_Si_Si = 0;
	repulsive_energy_O_O = 0;

	//initialize acceleration dynamic array to 0 
	for (int i = 0; i < totalAtomQuantity * 3; ++i)
		acceleration[i] = 0.0;


	//calculate length between each two bounded atoms (Si and O)
	int a = 1;
	for (int i = 0; i < atomQuantity[1]; ++i) 
	{
		//store Si position using atomCoordiante(need to know start point of Si atom)
		directCenterCoorSi[0] = coordinates[i * 3 + startPoint[1]];
		directCenterCoorSi[1] = coordinates[i * 3 + 1 + startPoint[1]];
		directCenterCoorSi[2] = coordinates[i * 3 + 2 + startPoint[1]];

		//store bonded O position 
		for (int j = 0; j < 4; ++j)
		{
			orderObondedwithSi[j] = groupSi[i][j];//get bonded O orders

			//get coordinates of each bonded O
			directBondedCoorOwithSi[j][0] = coordinates[orderObondedwithSi[j] * 3 + startPoint[0]];//x
			directBondedCoorOwithSi[j][1] = coordinates[orderObondedwithSi[j] * 3 + 1 + startPoint[0]];//y
			directBondedCoorOwithSi[j][2] = coordinates[orderObondedwithSi[j] * 3 + 2 + startPoint[0]];//z
		}

		//get real direct shortest coordinates of bonded O with center Si
		for (int j = 0; j < 4; ++j)
			EFTool.directShortestCoor(directCenterCoorSi, directBondedCoorOwithSi[j], realDirectBondedCoorOwithSi[j]);

			
		//--------------------------------------------------------------------------------------------------------------

		//convert Direct to Cartesian coordinates
		EFTool.directtoCart(directCenterCoorSi, cartesianCenterCoorSi);//convert center Si coordiantes from Direct to Cartesian


		for (int j = 0; j < 4; ++j)
			EFTool.directtoCart(realDirectBondedCoorOwithSi[j], realCartesianDirectBondedCoorOwithSi[j]);
		

		//length between center Si and bonded O
		T length_Si_O1 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0]);
		T length_Si_O2 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[1]);
		T length_Si_O3 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[2]);
		T length_Si_O4 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[3]);


		//angle
		T cosSi_O1_O2 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0], realCartesianDirectBondedCoorOwithSi[1], length_Si_O1, length_Si_O2);
		T cosSi_O1_O3 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0], realCartesianDirectBondedCoorOwithSi[2], length_Si_O1, length_Si_O3);
		T cosSi_O1_O4 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0], realCartesianDirectBondedCoorOwithSi[3], length_Si_O1, length_Si_O4);
		T cosSi_O2_O3 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[1], realCartesianDirectBondedCoorOwithSi[2], length_Si_O2, length_Si_O3);
		T cosSi_O2_O4 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[1], realCartesianDirectBondedCoorOwithSi[3], length_Si_O2, length_Si_O4);
		T cosSi_O3_O4 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[2], realCartesianDirectBondedCoorOwithSi[3], length_Si_O3, length_Si_O4);




		//--------------------------------------------------------------------------------------------------------------

		/*calculate radial energy and radial force
		 * and store forces to  */

		 //coefficient
		T diff_Si_O1 = length_Si_O1 - b0;
		T diff_Si_O2 = length_Si_O2 - b0;
		T diff_Si_O3 = length_Si_O3 - b0;
		T diff_Si_O4 = length_Si_O4 - b0;


		//radial energy
		radial_energy += 0.5 * kb *
			(diff_Si_O1*diff_Si_O1 + diff_Si_O2 * diff_Si_O2 +
				diff_Si_O3 * diff_Si_O3 + diff_Si_O4 * diff_Si_O4);

		//radial force for Si
		T fcoe = kb * diff_Si_O1 / length_Si_O1;
		Rf_O1_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[0][0] - cartesianCenterCoorSi[0]);//force on x by O1
		Rf_O1_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[0][1] - cartesianCenterCoorSi[1]);//force on y by O1
		Rf_O1_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[0][2] - cartesianCenterCoorSi[2]);//force on z by O1

		fcoe = kb * diff_Si_O2 / length_Si_O2;
		Rf_O2_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[1][0] - cartesianCenterCoorSi[0]);//force on x by O2
		Rf_O2_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[1][1] - cartesianCenterCoorSi[1]);//force on y by O2
		Rf_O2_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[1][2] - cartesianCenterCoorSi[2]);//force on z by O2

		fcoe = kb * diff_Si_O3 / length_Si_O3;
		Rf_O3_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[2][0] - cartesianCenterCoorSi[0]);//force on x by O3
		Rf_O3_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[2][1] - cartesianCenterCoorSi[1]);//force on y by O3
		Rf_O3_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[2][2] - cartesianCenterCoorSi[2]);//force on z by O3

		fcoe = kb * diff_Si_O4 / length_Si_O4;
		Rf_O4_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[3][0] - cartesianCenterCoorSi[0]);//force on x by O4
		Rf_O4_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[3][1] - cartesianCenterCoorSi[1]);//force on y by O4
		Rf_O4_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[3][2] - cartesianCenterCoorSi[2]);//force on z by O4

		//sum all the radial forces applied on Si
		Rf_on_Si[0] = Rf_O1_on_Si[0] + Rf_O2_on_Si[0] + Rf_O3_on_Si[0] + Rf_O4_on_Si[0];
		Rf_on_Si[1] = Rf_O1_on_Si[1] + Rf_O2_on_Si[1] + Rf_O3_on_Si[1] + Rf_O4_on_Si[1];
		Rf_on_Si[2] = Rf_O1_on_Si[2] + Rf_O2_on_Si[2] + Rf_O3_on_Si[2] + Rf_O4_on_Si[2];

		//store radial force into acceleration(need to know startpoint of Si)
		acceleration[i * 3 + startPoint[1]] += (Rf_on_Si[0] / massSi);
		acceleration[i * 3 + 1 + startPoint[1]] += (Rf_on_Si[1] / massSi);
		acceleration[i * 3 + 2 + startPoint[1]] += (Rf_on_Si[2] / massSi);

		//store radial force into acceleration which just add minus sign to Si forces(need to know startpoint of O)
		acceleration[orderObondedwithSi[0] * 3 + startPoint[0]] -= (Rf_O1_on_Si[0] / massO);//O1, x
		acceleration[orderObondedwithSi[0] * 3 + 1 + startPoint[0]] -= (Rf_O1_on_Si[1] / massO);//O1, y
		acceleration[orderObondedwithSi[0] * 3 + 2 + startPoint[0]] -= (Rf_O1_on_Si[2] / massO);//O1, z

		acceleration[orderObondedwithSi[1] * 3 + startPoint[0]] -= (Rf_O2_on_Si[0] / massO);//O2, x
		acceleration[orderObondedwithSi[1] * 3 + 1 + startPoint[0]] -= (Rf_O2_on_Si[1] / massO);//O2, y
		acceleration[orderObondedwithSi[1] * 3 + 2 + startPoint[0]] -= (Rf_O2_on_Si[2] / massO);//O2, z

		acceleration[orderObondedwithSi[2] * 3 + startPoint[0]] -= (Rf_O3_on_Si[0] / massO);//O3, x
		acceleration[orderObondedwithSi[2] * 3 + 1 + startPoint[0]] -= (Rf_O3_on_Si[1] / massO);//O3, y
		acceleration[orderObondedwithSi[2] * 3 + 2 + startPoint[0]] -= (Rf_O3_on_Si[2] / massO);//O3, z

		acceleration[orderObondedwithSi[3] * 3 + startPoint[0]] -= (Rf_O4_on_Si[0] / massO);//O4, x
		acceleration[orderObondedwithSi[3] * 3 + 1 + startPoint[0]] -= (Rf_O4_on_Si[1] / massO);//O4, y
		acceleration[orderObondedwithSi[3] * 3 + 2 + startPoint[0]] -= (Rf_O4_on_Si[2] / massO);//O4, z



		//--------------------------------------------------------------------------------------------------------------

		/*calculate angular energy and angular force
		 * and store forces to acceleration */

		//coefficient
		T diff_Si_O1_O2 = cosSi_O1_O2 - cos0_Si;
		T diff_Si_O1_O3 = cosSi_O1_O3 - cos0_Si;
		T diff_Si_O1_O4 = cosSi_O1_O4 - cos0_Si;
		T diff_Si_O2_O3 = cosSi_O2_O3 - cos0_Si;
		T diff_Si_O2_O4 = cosSi_O2_O4 - cos0_Si;
		T diff_Si_O3_O4 = cosSi_O3_O4 - cos0_Si;

		//angular energy for Si
		angular_energy_Si += 0.5 * kthetaSi * b0 * b0 *
			(diff_Si_O1_O2*diff_Si_O1_O2 + diff_Si_O1_O3 * diff_Si_O1_O3 +
				diff_Si_O1_O4 * diff_Si_O1_O4 + diff_Si_O2_O3 * diff_Si_O2_O3 +
				diff_Si_O2_O4 * diff_Si_O2_O4 + diff_Si_O3_O4 * diff_Si_O3_O4);

		//angular force for Si
		fcoe = -1.0 * kthetaSi * b0 * b0;

		/* force calculation (Si, O1 and O2) */
		for (int m = 0; m < 3; m++)//force between Si O1 and O2
			fSiO1O2_on_O1[m] = fcoe * ((diff_Si_O1_O2 / length_Si_O1) *
			((realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2 
				- cosSi_O1_O2 * (realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1));
		for (int m = 0; m < 3; m++)//force between Si O1 and O2
			fSiO1O2_on_O2[m] = fcoe * ((diff_Si_O1_O2 / length_Si_O2) *
			((realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1 
				- cosSi_O1_O2 * (realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
			fSiO1O2_on_Si[m] = -(fSiO1O2_on_O1[m] + fSiO1O2_on_O2[m]);


		/* force calculation (Si, O1 and O3) */
		for (int m = 0; m < 3; m++)//force between Si O1 and O3
			fSiO1O3_on_O1[m] = fcoe * ((diff_Si_O1_O3 / length_Si_O1) *
			((realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3 
				- cosSi_O1_O3 * (realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1));
		for (int m = 0; m < 3; m++)//force between Si O1 and O3
			fSiO1O3_on_O3[m] = fcoe * ((diff_Si_O1_O3 / length_Si_O3) *
			((realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1 
				- cosSi_O1_O3 * (realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O3 on Si
			fSiO1O3_on_Si[m] = -(fSiO1O3_on_O1[m] + fSiO1O3_on_O3[m]);


		/* force calculation (Si, O1 and O4) */
		for (int m = 0; m < 3; m++)//force between Si O1 and O4
			fSiO1O4_on_O1[m] = fcoe * ((diff_Si_O1_O4 / length_Si_O1) *
			((realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4 
				- cosSi_O1_O4 * (realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1));
		for (int m = 0; m < 3; m++)//force between Si O1 and O4
			fSiO1O4_on_O4[m] = fcoe * ((diff_Si_O1_O4 / length_Si_O4) *
			((realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1 
				- cosSi_O1_O4 * (realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
			fSiO1O4_on_Si[m] = -(fSiO1O4_on_O1[m] + fSiO1O4_on_O4[m]);


		/* force calculation (Si, O2 and O3) */
		for (int m = 0; m < 3; m++)//force between Si O2 and O3
			fSiO2O3_on_O2[m] = fcoe * ((diff_Si_O2_O3 / length_Si_O2) *
			((realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3 
				- cosSi_O2_O3 * (realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2));
		for (int m = 0; m < 3; m++)//force between Si O2 and O3
			fSiO2O3_on_O3[m] = fcoe * ((diff_Si_O2_O3 / length_Si_O3) *
			((realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2 
				- cosSi_O2_O3 * (realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O3 on Si
			fSiO2O3_on_Si[m] = -(fSiO2O3_on_O2[m] + fSiO2O3_on_O3[m]);


		/* force calculation (Si, O2 and O4) */
		for (int m = 0; m < 3; m++)//force between Si O2 and O4
			fSiO2O4_on_O2[m] = fcoe * ((diff_Si_O2_O4 / length_Si_O2) *
			((realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4 
				- cosSi_O2_O4 * (realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2));
		for (int m = 0; m < 3; m++)//force between Si O2 and O4
			fSiO2O4_on_O4[m] = fcoe * ((diff_Si_O2_O4 / length_Si_O4) *
			((realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2 
				- cosSi_O2_O4 * (realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O4 on Si
			fSiO2O4_on_Si[m] = -(fSiO2O4_on_O2[m] + fSiO2O4_on_O4[m]);


		/* force calculation (Si, O3 and O4) */
		for (int m = 0; m < 3; m++)//force between Si O3 and O4
			fSiO3O4_on_O3[m] = fcoe * ((diff_Si_O3_O4 / length_Si_O3) *
			((realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4 
				- cosSi_O3_O4 * (realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3));
		for (int m = 0; m < 3; m++)//force between Si O3 and O4
			fSiO3O4_on_O4[m] = fcoe * ((diff_Si_O3_O4 / length_Si_O4) *
			((realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3 
				- cosSi_O3_O4 * (realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O3 and O4 on Si
			fSiO3O4_on_Si[m] = -(fSiO3O4_on_O3[m] + fSiO3O4_on_O4[m]);

		//store angular force into forceAcceleration(need to know startpoint of Si)
		for (int m = 0; m < 3; m++)
			acceleration[i * 3 + m + startPoint[1]] += (fSiO1O2_on_Si[m] + fSiO1O3_on_Si[m] + fSiO1O4_on_Si[m] +
				fSiO2O3_on_Si[m] + fSiO2O4_on_Si[m] + fSiO3O4_on_Si[m]) / massSi;


		//store angular force into forceAcceleration(need to know startpoint of O)
		for (int m = 0; m < 3; m++) {
			acceleration[orderObondedwithSi[0] * 3 + m + startPoint[0]] += ((fSiO1O2_on_O1[m] + fSiO1O3_on_O1[m] + fSiO1O4_on_O1[m]) / massO);
			acceleration[orderObondedwithSi[1] * 3 + m + startPoint[0]] += ((fSiO1O2_on_O2[m] + fSiO2O3_on_O2[m] + fSiO2O4_on_O2[m]) / massO);
			acceleration[orderObondedwithSi[2] * 3 + m + startPoint[0]] += ((fSiO1O3_on_O3[m] + fSiO2O3_on_O3[m] + fSiO3O4_on_O3[m]) / massO);
			acceleration[orderObondedwithSi[3] * 3 + m + startPoint[0]] += ((fSiO1O4_on_O4[m] + fSiO2O4_on_O4[m] + fSiO3O4_on_O4[m]) / massO);
		}


		// calculate repulsive energy and force between Si and Si
		for (int k = a; k < atomQuantity[1]; ++k)
		{
			//store all other Si coordiante in Direct(need to know startpoint of Si)
			directRepSiSi[0] = coordinates[k * 3 + startPoint[1]];
			directRepSiSi[1] = coordinates[k * 3 + 1 + startPoint[1]];
			directRepSiSi[2] = coordinates[k * 3 + 2 + startPoint[1]];

			//get real direct coor of Si near the center Si
			EFTool.directShortestCoor(directCenterCoorSi, directRepSiSi, directRepSiSishortest);

			//convert direct to cartesian 
			EFTool.directtoCart(directRepSiSishortest, cartesianRepSiSishortest);

			//get length between two Si atoms using Cartesian	
			repSiSilength = EFTool.distanceinCart(cartesianCenterCoorSi, cartesianRepSiSishortest);

			if (repSiSilength <= 3.10)
			{
				//calculate energy
				repulsive_energy_Si_Si += 0.5 * repulsiveCoeforSiSi * (3.10 - repSiSilength)*(3.10 - repSiSilength)*(3.10 - repSiSilength);

				//calculate force
				for (int m = 0; m < 3; ++m)
				{
					frepSiSi[m] = -3 * repulsiveCoeforSiSi *
						(3.10 - repSiSilength)*(3.10 - repSiSilength) / repSiSilength *
						(cartesianRepSiSishortest[m] - cartesianCenterCoorSi[m]);
					//store force in acceleration (need to know startpoint of Si)
					acceleration[i * 3 + m + startPoint[1]] += (frepSiSi[m] / massSi);
					acceleration[k * 3 + m + startPoint[1]] -= (frepSiSi[m] / massSi);
				}
			}
		}
		++a;


	}


	//--------------------------------------------------------------------------------------------------------------



	//calculate angular energy of O
	int b = 1;
	for (int i = 0; i <atomQuantity[0]; ++i)
	{
		//store O position(need to know the startpoint of O)
		directCenterCoorO[0] = coordinates[i * 3 + startPoint[0]];
		directCenterCoorO[1] = coordinates[i * 3 + 1 + startPoint[0]];
		directCenterCoorO[2] = coordinates[i * 3 + 2 + startPoint[0]];

		//get order of Si bonded with center O
		for (int j = 0; j < 2; ++j)
		{
			orderSibondedwithO[j] = groupO[i][j];
			for (int k = 0; k < 3; ++k)//get direct coordinates of Si(need to know the startpoint of SI)
				directBondedCoorSiwithO[j][k] = coordinates[orderSibondedwithO[j] * 3 + k + startPoint[1]];
		}

		//get real direct coordinate of Si that bonded with O
		for (int j = 0; j < 2; ++j)
			EFTool.directShortestCoor(directCenterCoorO, directBondedCoorSiwithO[j], realDirectBondedCoorSiwithO[j]);

		//convert coordinates from Direct to Cartesian
		EFTool.directtoCart(directCenterCoorO, cartesianCenterCoorO);//for center O
		for (int j = 0; j < 2; ++j)//for bonded Si
			EFTool.directtoCart(realDirectBondedCoorSiwithO[j], realCartesianDirectBondedCoorSiwithO[j]);

		//length
		T length_O_Si1 = EFTool.distanceinCart(cartesianCenterCoorO, realCartesianDirectBondedCoorSiwithO[0]);
		T length_O_Si2 = EFTool.distanceinCart(cartesianCenterCoorO, realCartesianDirectBondedCoorSiwithO[1]);

		//angle
		T cosO_Si1_Si2 = EFTool.angleinCart(cartesianCenterCoorO, realCartesianDirectBondedCoorSiwithO[0], realCartesianDirectBondedCoorSiwithO[1], length_O_Si1, length_O_Si2);

		//calculate angular energy
		T diff_O_Si1_Si2 = cosO_Si1_Si2 - cos0_O;
		angular_energy_O += 0.5 * kthetaO * b0 * b0 * diff_O_Si1_Si2 * diff_O_Si1_Si2;

		//calculate angular force
		T fcoe = -1.0 * kthetaO * b0 * b0;
		//force calculation (Si, O2 and O4)
		for (int m = 0; m < 3; m++)//force between O Si1 and Si2
			fOSi1Si2_on_Si1[m] = fcoe * ((diff_O_Si1_Si2 / length_O_Si1) *
			((realCartesianDirectBondedCoorSiwithO[1][m] - cartesianCenterCoorO[m]) / length_O_Si2 
				- cosO_Si1_Si2 * (realCartesianDirectBondedCoorSiwithO[0][m] - cartesianCenterCoorO[m]) / length_O_Si1));
		for (int m = 0; m < 3; m++)//force between Si O2 and O4
			fOSi1Si2_on_Si2[m] = fcoe * ((diff_O_Si1_Si2 / length_O_Si2) *
			((realCartesianDirectBondedCoorSiwithO[0][m] - cartesianCenterCoorO[m]) / length_O_Si1 
				- cosO_Si1_Si2 * (realCartesianDirectBondedCoorSiwithO[1][m] - cartesianCenterCoorO[m]) / length_O_Si2));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O4 on Si
			fOSi1Si2_on_O[m] = -(fOSi1Si2_on_Si1[m] + fOSi1Si2_on_Si2[m]);

		//store angular force into atomAccelelation(need to know the startpoint of O)
		for (int m = 0; m < 3; m++)
			acceleration[i * 3 + m + startPoint[0]] += (fOSi1Si2_on_O[m] / massO);

		//store angular force into atomAccelelation(need to know the startpoint of Si)
		for (int m = 0; m < 3; m++)
		{
			acceleration[orderSibondedwithO[0] * 3 + m + startPoint[1]] += (fOSi1Si2_on_Si1[m] / massSi);
			acceleration[orderSibondedwithO[1] * 3 + m + startPoint[1]] += (fOSi1Si2_on_Si2[m] / massSi);
		}



		//calculate repulsive energy between O and O
		for (int j = b; j < atomQuantity[0]; ++j)//get direct coordinates of center O
		{
			directRepOO[0] = coordinates[j * 3 + startPoint[0]];
			directRepOO[1] = coordinates[j * 3 + 1 + startPoint[0]];
			directRepOO[2] = coordinates[j * 3 + 2 + startPoint[0]];

			//find real direct coordinates of O near center O
			EFTool.directShortestCoor(directCenterCoorO, directRepOO, directRepOOshortest);

			//convert coordiantes from Direct to Cartesian for RepOO
			EFTool.directtoCart(directRepOOshortest, cartesianRepOOshortest);

			//get length
			repOOlength = EFTool.distanceinCart(cartesianCenterCoorO, cartesianRepOOshortest);

			if (repOOlength <= 2.533)
			{
				//calculate energy
				repulsive_energy_O_O += 0.5*repulsiveCoeforOO*(2.533 - repOOlength)*(2.533 - repOOlength)*(2.533 - repOOlength);

				//calculate force for O and O
				for (int m = 0; m < 3; ++m)
				{
					frepOO[m] = -3 * repulsiveCoeforOO * (2.533 - repOOlength)*(2.533 - repOOlength) / repOOlength *
						(cartesianRepOOshortest[m] - cartesianCenterCoorO[m]);
					//store acceleration to acceleration(need to know the startpoint of O)
					acceleration[i * 3 + m + startPoint[0]] += (frepOO[m] / massO);
					acceleration[j * 3 + m + startPoint[0]] -= (frepOO[m] / massO);
				}
			}
		}
		++b;

	}

	//--------------------------------------------------------------------------------------------------------------

	/*std::cout << "radial energy: " << radial_energy << std::endl;
	std::cout << "Si angular energy: " << angular_energy_Si << std::endl;
	std::cout << "O angular energy: " << angular_energy_O << std::endl;
	std::cout << "Si repulsive energy: " << repulsive_energy_Si_Si << std::endl;
	std::cout << "O repulsive energy: " << repulsive_energy_O_O << std::endl;*/

	total_energy = radial_energy + angular_energy_Si + angular_energy_O +
		repulsive_energy_O_O + repulsive_energy_Si_Si;



};
