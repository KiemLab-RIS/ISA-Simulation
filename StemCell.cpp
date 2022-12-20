//
//  StemCell.cpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#include "StemCell.hpp"
//
//
//
StemCell::StemCell(int cloneID,int age,double ef,int recTime) {
	_cellID = CellID++;
	_cloneID = cloneID;
	_alive = 1;
	_cycleCount = 0;
	_age = age;
	_state = 99;
	_epigeneticFactor = ef;
	_recoveryTime = recTime;
	_pStem=0.0;
	_pDif=0.0;
	_pApop=0.0;
	_pCycle=0.0;
}
