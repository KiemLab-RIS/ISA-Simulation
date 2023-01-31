//
//  main.cpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//
//
//  This program simulates bone marrow HSC behavior and
//  also does a very simplified simulation of blood production
//  from HSCs. This simulation expands blood cells
//  of various types (gran/mono/nk etc), with a bit of time delay, into a
//  collections (very) roughly cooresponding to a blood cell population that
//  we can then simulate performing Insertion Site Analysis (ISA) processing on.
//
//  This simulation is only meant to show how time delay expansion of progenitor
//  cells can be sampled with ISA. There is no claim that any of these expansion
//  pathways are actually realistic. We don't have that information in enough detail
//  yet that it would be worthwhile to simulate
//
//  To run this program, make an entry in the TEST_PARAM array
//  -------------------------------------------------------------------
//	int cells;      starting cells
//	double r1;      initial rep rate
//	double r2;      long-term rep rate
//	int day;        day to switch from initial to long term....if used
//	double pStem;   chance HSC daughter stays stem
//	double epiMean; HSC cycle time mean
//	double epiVar;  HSC cycle time var
//	std::string s; store results
//
//  for example
//
//	{120000,0.8,0.8,20,0.52,1.0,0.05,"/Volumes/Fred Hutch Data/LineageSimResults/Z_120000_8_8_20_6_1_5/"},
//	120000  - starting stem cells
//	0.8     - HSC replication rate at start (is really inverse...20% cycle per 'day')
//	0.8     - HSC replication rate for long term...since both 0.8, 20% cycle...not change
//	20      - day cell cycle would change
//	0.52    - chance stay stem when hsc division occurs
//	1.0     - mean for epigenetic cycle time differences
//	0.05    - var for epigenetic cycle differnce
//	s       - store results
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
