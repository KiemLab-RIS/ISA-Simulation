//
//  BloodCell.hpp
//  StemSim2
//
//  Created by mark enstrom on 4/14/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#ifndef BloodCell_hpp
#define BloodCell_hpp

#include <stdio.h>
#include <vector>
#include <random>
#include <iostream>
#include "DiffCell.hpp"

class CloneDev {
public:
	int _cloneID;
	int _count;
	int _type;
};



class BloodCell {
public:
	int _cloneID;
	int	_cellID;	
	int _count;
	int _age;
	int _type;
	
	BloodCell(int cloneID,int cellID,int count,int type);
};

typedef std::vector<CloneDev*> VFIFO;

#define BLOOD_FIFO_LENGTH 1000
#define BLOOD_UNIT_SIZE 100000

class BloodFIFO {
public:
	int _head;
	int _length = BLOOD_FIFO_LENGTH;
	VFIFO *_p;
	BloodFIFO() {
		_p = new VFIFO[BLOOD_FIFO_LENGTH];
	}
	//
	// expand each clone cloneOutput times and:
	//     	break the clone cells into equal units with 50k? cells each
	//		schedule the cells to mature using a gaussion distribution
	//      based on delay(days) and var(variance)
	//
	void distribute(std::vector<DiffCell*> *pvCells,
					int fixedDelay,
					int expansionFactor,
					double eVar,  //expansion variance
					double fifoVar) // maturation variance
	{
		long cellsAdded = 0;

		for (auto pCell:*pvCells) {
			//
			// epigenetic factor for cell =
			// normal distribution(0.5,0.04);
			//
			double tFactor = expansionFactor * pCell->_ef;
			if (tFactor > 22.0) {
				tFactor = 22.0;
			}
			int cloneOutput = std::pow(2,tFactor);
			
			//std::cout << "exp = " << expansionFactor << "   ef = " << pCell->_ef << " clone output = " << cloneOutput << "\n";
			
			cellsAdded += cloneOutput;
			
			int i = 1;
			while (cloneOutput > 0) {
				int delay = fixedDelay+i;
				int number = cloneOutput;
				if (number > BLOOD_UNIT_SIZE) {
					number = BLOOD_UNIT_SIZE;
				}
				// cells mature after delay
				this->add(pCell->_cloneID,number,pCell->_type,delay);
				//std::cout << "add cell ID [" << pCell->_cloneID << "] rf = " << rf << " cells =  " << number << " at delay = " << delay << std::endl;
				i++;
				cloneOutput = cloneOutput - number;
			}
		}
		for (int i  = 0; i < pvCells->size(); i++) {
			delete (*pvCells)[i];
		}
		pvCells->clear();
		std::cout << "      Cells added after expansion = " << cellsAdded << "\n";
	}
	//
	// add cells to FIFO at given day delay
	//
	void add(int cloneID,int count,int type,int day) {
		if (cloneID > 1200000) {
			std::cout << "bad clone ID " << cloneID << std::endl;
		}
		int insertDay = _head + day;
		if (insertDay >= _length) {
			insertDay = insertDay - _length;
		}
		CloneDev *cd = new CloneDev();
		cd->_cloneID = cloneID;
		cd->_count = count;
		cd->_type = type;
		_p[insertDay].push_back(cd);
	}
	//
	// get all cells that mature at the head of the FIFO, clear and advance FIFO
	//
	void get(std::vector<BloodCell*> *bloodCells) {
		for (CloneDev *pr:_p[_head]) {
			if (pr->_cloneID > 1200000) {
				std::cout << "ERROR bad clone ID \n";
			}
			BloodCell *pGranBlood = new BloodCell(pr->_cloneID,88,pr->_count,pr->_type);
			bloodCells->push_back(pGranBlood);
			delete pr;
		}
		_p[_head].clear();
		_head ++;
		if (_head >= _length) {
			_head = 0;
		}
		return;
	}
	
	
	void Diag() {
		std::cout << "head = " << _head << std::endl;
		for (int i = 0; i < _length; i++) {
			auto v = _p[i];
			std::cout << "size at " << i << " = " << v.size() << std::endl;
		}
	}
};


#endif /* BloodCell_hpp */
