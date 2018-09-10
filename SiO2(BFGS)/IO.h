#ifndef IO_H_
#define IO_H_

#include <fstream>
#include <iostream>
#include "CONSTANT.h"



class IO
{
public:
	IO();
	~IO();

public:
	void readfromFile();
	void output(const int &n);
};

#endif /* IO_H_ */
