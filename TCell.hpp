//
//  TCell.hpp
//  StemMPP
//
//  Created by mark enstrom on 4/30/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef TCell_hpp
#define TCell_hpp

#include <stdio.h>
#include "DiffCell.hpp"


class TCell:public DiffCell {
public:
	TCell(int cloneID,int cellID,int age,double ef)
		: DiffCell(cloneID,cellID,1,0,0,age,T_TYPE,ef)
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

#endif /* TCell_hpp */
