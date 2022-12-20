//
//  GranCell.hpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef GranCell_hpp
#define GranCell_hpp

#include <stdio.h>
#include "DiffCell.hpp"

class GranCell:public DiffCell {
public:
	GranCell(int cloneID,int cellID,int age,double ef)
		: DiffCell(cloneID,cellID,1,0,0,age,GRAN_TYPE,ef)
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




#endif /* GranCell_hpp */
