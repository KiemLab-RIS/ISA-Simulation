//
//  MppCell.hpp
//  StemMPP
//
//  Created by mark enstrom on 4/23/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef MppCell_hpp
#define MppCell_hpp

#include <stdio.h>

enum MppState {rest,cycle,diff};

class MppCell {
public:
	int _cloneID;
	float _repTime;
	MppState _state;
	int _cycleCount;
	int _alive;
	int _age;
	double _ef;
	MppCell(int cloneID,int age,double ef);
};
#endif /* MppCell_hpp */
