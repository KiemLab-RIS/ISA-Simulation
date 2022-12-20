//
//  MppCell.cpp
//  StemMPP
//
//  Created by mark enstrom on 4/23/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#include "MppCell.hpp"

MppCell::MppCell(int cloneID,int age,double ef) {
	_cloneID = cloneID;
	_alive = 1;
	_cycleCount = 0;
	_age = age;
	_state = rest;
	_ef = ef;
}
