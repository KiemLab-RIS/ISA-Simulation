//
//  main.cpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <iterator>
#include <algorithm>
#include <list>
#include <map>
#include <deque>
#include <fstream>
#include <chrono>
#include <set>
#include <filesystem>
#include "unordered_set"
//#include "StemCell.hpp"
//#include "GranCell.hpp"
//#include "BloodCell.hpp"
//#include "MonoCell.hpp"
//#include "BCell.hpp"
//#include "TCell.hpp"
//#include "NkCell.hpp"
//#include "MppCell.hpp"

#include "Simulation.hpp"


int StemCell::CellID = 0;
using namespace std;

#define STEM_CELL_CYCLE_TIME 1



int main(int argc, const char * argv[]) {
	
	//   Starting Cells
	//   initial replication rate
	//   long-term replication rate
	//   short-long switch
	//   pStem daughter cell 1
	//   pStem daughter cell 2
	//   death rate
	//   dst dir
	//
	TEST_PARAM tp[] = {

		//{120000,0.0,0.45,20,0.5,0.5,0.3,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.43,20,0.5,0.5,0.3,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.43,10,0.5,0.5,0.3,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.43,10,0.5,0.5,0.2,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.43,10,0.5,0.5,0.15,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.41,10,0.5,0.5,0.15,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.41,10,0.5,0.5,0.25,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.41,20,0.5,0.5,0.25,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.41,20,0.5,0.5,0.20,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.41,20,0.5,0.5,0.15,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.0,0.41,20,0.5,0.5,0.10,"/Volumes/Fred Hutch Data/LineageSimResults/Z_80000_cp_1/"},
		//{120000,0.9,0.9,20,0.5,1.0,0.00,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_9_9_20_1_0/"},
		//{120000,0.9,0.9,20,0.5,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_9_9_20_1_5/"},
		//{120000,0.8,0.8,20,0.5,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_8_8_20_1_5/"},
		//{120000,0.0,0.8,20,0.5,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_0_8_20_1_5/"},
		//{120000,0.0,0.9,20,0.5,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_0_9_20_1_05/"},

		//{120000,0.0,0.9,20,0.5,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_0_9_20_1_05/"},
		//
		// a large clone developed!
		//
		{120000,0.8,0.8,20,0.6,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_8_8_20_6_1_5/"},
		//
		//
		//
		{120000,0.8,0.8,20,0.52,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_8_8_20_52_1_5/"},
		
		
		
	};
	
	for (auto test:tp) {
		Simulation sim = Simulation(&test);
		sim.run_sim();
	}
	

}
