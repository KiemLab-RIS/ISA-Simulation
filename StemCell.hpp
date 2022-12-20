//
//  StemCell.hpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef StemCell_hpp
#define StemCell_hpp

#include <stdio.h>
//5



//
// 
//
//
//
class StemCell {
public:
	static int CellID;
	int _cloneID;
	int	_cellID;
	int _state;
	int _cycleCount;
	int _recoveryCount;
	int _recoveryTime;
	double _epigeneticFactor;
	int _alive;
	int _age;
	
	double _pStem;
	double _pDif;
	double _pApop;
	double _pCycle;
	
	StemCell(int cloneID,int age,double ef,int recTim);
};


#endif /* StemCell_hpp */
