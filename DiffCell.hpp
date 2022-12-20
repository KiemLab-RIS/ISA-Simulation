//
//  DiffCell.hpp
//  StemSim2
//
//  Created by mark enstrom on 4/16/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef DiffCell_hpp
#define DiffCell_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>

#define GRAN_TYPE 1
#define MONO_TYPE 2
#define B_TYPE 3
#define T_TYPE 4
#define NK_TYPE 5


class DiffCell {
public:
	int _cloneID;
	int	_cellID;
	int _count;
	int _state;
	int _maturity;
	int _age;
	int _type;
	double _ef;
	DiffCell(int cloneid, int cellid,int count,int state,int mat,int age,int type,double ef)
	:	_cloneID{cloneid},
	_cellID{cellid},
	_count{count},
	_state{state},
	_maturity{mat},
	_age{age},
	_type{type},
	_ef{ef}
	{}
public:
	void advanceMaturity() {_maturity++;}

};
//
// helper function
//

void doSample(std::vector<int> cellUnits,
			  std::string cellType,
			  std::string sDate,
              std::string outputDirectory);

void doSampleAll(std::map<int,double> &cloneCounts,
				 std::string cellType,
				 std::string sDate,
				 std::string outputDirectory);


#endif /* DiffCell_hpp */
