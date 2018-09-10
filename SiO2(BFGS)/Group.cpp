#include "Group.h"



Group::Group()
{
	//std::cout << "Group constructor" << std::endl;
}


Group::~Group()
{
	//std::cout << "Group deconstructor" << std::endl;
}


/*find all bonded atoms*/
void Group::Graph()
{
	Tool groupTool;

	//find place of atom O in POSCAR
	for (int i = 0; i < atomType; ++i)
	{
		if (atomName[i] == "O")
		{
			orderofO = i;//get the order of O in POSCAR
			break;
		}
	}

	//get each real O atom coordinate x, y, x from atomCoordinate
	for (int i = startPoint[orderofO]; i < atomQuantity[orderofO] * 3; i += 3)// i is NOT the order of O
	{	
		directCoorO[0] = atomCoordinate[i];//get x
		directCoorO[1] = atomCoordinate[i + 1];//get y
		directCoorO[2] = atomCoordinate[i + 2];//get z

		//get all possible coordinates of this O atom
		allPossibleCoordinates(directCoorO, allCoorO);

		//calcualte distance between Si and O, P and O;
		for (int j = startPoint[1]; j < totalAtomQuantity * 3; j += 3)//get Si coordinates(j is NOT equal to the order of Si)
		{
			directCoorSi[0] = atomCoordinate[j];//get x
			directCoorSi[1] = atomCoordinate[j + 1];//get y
			directCoorSi[2] = atomCoordinate[j + 2];//get z

			//convert direct to cartesian for both Si and O coordinates
			groupTool.directtoCart(directCoorSi, cartesianCoorSi);//for Si

			for (int k = 0; k < 27; ++k)//find shortest distance between Si and 27 possible O 
			{
				
				directCoorO[0] = allCoorO[k][0];//get O coordinate of x
				directCoorO[1] = allCoorO[k][1];//get O coordinate of y
				directCoorO[2] = allCoorO[k][2];//get O coordinate of z
				groupTool.directtoCart(directCoorO, cartesianCoorO);//for O

				T distance = groupTool.distanceinCart(cartesianCoorO, cartesianCoorSi);//get distance

				//determine whether they bond together
				int temp = 0;
				if (distance < 1.56)//the bond Slength between Si and O is around 1.55A
				{
					//[(j - startPoint[1]) / 3] is the order of Si, [i / 3] is the order of O.
					graph[(j - startPoint[1]) / 3][i / 3] = 1;//if Si and O bonded together, set value = 1.
				}

			}
			

		}
	}
}


/*get all possible coordinates for one atom*/
void Group::allPossibleCoordinates(const T (&coor)[3], T (&allCoor)[27][3])
{
	int loop = 0;
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			for (int k = -1; k <= 1; k++)
			{
				allCoor[loop][0] = coor[0] + i;//x
				allCoor[loop][1] = coor[1] + j;//y
				allCoor[loop][2] = coor[2] + k;//z
				++loop;
			}
		}
	}
}


//search bonded group and store orders in each array
void Group::groupfromGraph()
{
	//group Si and store O's order in array
	for (int i = 0; i < atomQuantity[1]; ++i)
	{
		int n = 0;
		for (int j = 0; j < atomQuantity[0]; ++j) {
			if (graph[i][j] != 0 && n < 4) {
				//store O's order
				groupSi[i][n] = j;
				++n;
			}
		}
	}

	//group O and store Si's order in array
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		int m = 0;
		for (int j = 0; j < atomQuantity[1]; ++j)
		{
			if (graph[j][i] != 0 && m < 2)
			{
				groupO[i][m] = j;
				++m;
			}
		}
	}
};


//cut bonds
void Group::cutBond(int graph[][128])
{
	//choose a O atom as center
	int center_order = generator_int(0, atomQuantity[0] - 1);//need to know the quantity of O

	//find two Si atoms bounded to center O
	int a = 0;
	for (int i = 0; i < atomQuantity[1]; ++i)//need to know the quantity of Si
	{
		if (graph[i][center_order] != 0 && a < 2)
		{
			orderSi[a] = i;
			++a;
		}
	}
	//get other three O atoms bounded to two Si atoms
	int b1 = 0;
	int b2 = 0;
	for (int i = 0; i < atomQuantity[0]; ++i)//need to know the quantity of O
	{
		if (i != center_order)
		{
			if (graph[orderSi[0]][i] != 0 && b1 < 3)
			{
				ordernotcenterO[0][b1] = i;
				++b1;
			}
			if (graph[orderSi[1]][i] != 0 && b2 < 3)
			{
				ordernotcenterO[1][b2] = i;
				++b2;
			}

		}
	}

	//choose two O atoms randomly
	//for first Si
	int random1 = generator_int(0, 2);
	int orderO1 = ordernotcenterO[0][random1];
	//for second Si
	int random2 = generator_int(0, 2);
	int orderO2 = ordernotcenterO[1][random2];
	//exchange bonds
	graph[orderSi[0]][orderO1] = 0;
	graph[orderSi[1]][orderO2] = 0;

	graph[orderSi[0]][orderO2] = 1;
	graph[orderSi[1]][orderO1] = 1;


	//cout << "center O: "<<center_order+1<<endl;
	//cout << "Si1 and O1: "<< orderSi[0]+1<< "    "<<orderO1+1<<endl;
	//cout << "Si2 and O2: "<< orderSi[1]+1<< "    "<<orderO2+1<<endl;

};



//generate a integer in range (begin, end)
int Group::generator_int(int begin, int end)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(begin, end);
	return dis(gen);
}



