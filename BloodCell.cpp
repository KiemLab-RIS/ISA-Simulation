//
//  BloodCell.cpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#include "BloodCell.hpp"
BloodCell::BloodCell(int cloneID,int cellID,int count,int type) {
	_cellID = cellID;
	_cloneID = cloneID;
	_age = 0;
	_count = count;
	_type = type;
}
