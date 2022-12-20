//
//  MonoCell.hpp
//  StemSim2
//
//  Created by mark enstrom on 4/15/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef MonoCell_hpp
#define MonoCell_hpp

#include <stdio.h>
#include "DiffCell.hpp"


class MonoCell:public DiffCell {
public:
	MonoCell(int cloneID,int cellID,int age,double ef)
		: DiffCell(cloneID,cellID,1,0,0,age,MONO_TYPE,ef)
	{
		
	}
	int cloneID() {return _cloneID;}
	int age() {return _age;}
	int type() {return _type;}
	int count() {return _count;}
	void count(int c) {_count = c;}
	int maturity() {return _maturity;}
	void maturity(int m) {_maturity = m;}

};


#endif /* MonoCell_hpp */
