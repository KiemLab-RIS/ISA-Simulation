//
//  Simulation.hpp
//  StemCC
//
//  Created by mark enstrom on 6/29/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef Simulation_hpp
#define Simulation_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include "Simulation.hpp"
#include "StemCell.hpp"
#include "GranCell.hpp"
#include "BloodCell.hpp"
#include "MonoCell.hpp"
#include "BCell.hpp"
#include "TCell.hpp"
#include "NkCell.hpp"
#include "MppCell.hpp"

#define STEM_CELL_CYCLE_TIME 1

typedef struct _TEST_PARAM {
	int cells;
	double r1;
	double r2;
	int day;
	double pStem;
	double epiMean;
	double epiVar;
	std::string s;
}TEST_PARAM;


class Simulation {
public:
	TEST_PARAM* _tp;
	int day;
	int StartingStemCells;
	int _apopLimit = 600000;
	std::vector<std::pair<int,int>> uniqueCloneCount;
	std::vector<StemCell*> stemCellsResting;
	std::vector<StemCell*> stemCellsRecovering;
	std::vector<StemCell*> stemCellsCycling;
	std::vector<DiffCell*> granCells;
	std::vector<DiffCell*> monoCells;
	std::vector<DiffCell*> bCells;
	std::vector<DiffCell*> tCells;
	std::vector<DiffCell*> nkCells;
	std::vector<MppCell*>  mppCells;
	std::vector<BloodCell*> bloodCells;
	
	BloodFIFO *bfGran;
	BloodFIFO *bfMono;
	BloodFIFO *bfB;
	BloodFIFO *bfT;
	BloodFIFO *bfNK;

	
	Simulation(TEST_PARAM *ptp);
	int doSampleStem();
	void doSampleBlood();
	void simStemCellsV1();
	void simMppCellsV1();
	void simBloodCells();
	void run_sim();
};

#endif /* Simulation_hpp */
