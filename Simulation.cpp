//
//  Simulation.cpp
//  StemCC
//
//  Created by mark enstrom on 6/29/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//

#include <stdio.h>
#include "Simulation.hpp"

using namespace std;
//************************************************************************************
// map sort construct
//
//
//
//
//
//************************************************************************************
struct myclassFred {
    bool operator() (std::pair<int,int> pa1,std::pair<int,int> pa2) { return (pa1.second > pa2.second);}
} fredCompObj;



//************************************************************************************
//
//
//
//
//
//
//
//************************************************************************************
int Simulation::doSampleStem()
{
	char buffer[100];
	sprintf(buffer,"%04i",day);
	string sDate(buffer);
	//
	// summarize stem cells
	//
	std::vector<int> stemID(_tp->cells,0);
	for (StemCell* pCell : stemCellsResting) {
		stemID[pCell->_cloneID] += 1;
	}
	//
	// add cycling
	//
	for (StemCell* pCell : stemCellsCycling) {
		stemID[pCell->_cloneID] += 1;
	}
	//
	// add recovering
	//
	for (StemCell* pCell : stemCellsRecovering) {
		stemID[pCell->_cloneID] += 1;
	}
	
	std::vector<std::pair<int,int>> st;
	//
	// sort
	//
	for (int i = 0; i < _tp->cells; i++) {
		int nCells = stemID[i];
		if (nCells > 0) {
			st.push_back(std::pair<int,int>(i,nCells));
		}
	}
	std::sort(st.begin(),st.end(),fredCompObj);
	std::ofstream myfile;
	std::string filename = _tp->s +  "xSim_" + sDate;
	myfile.open (filename + "_stem_actual.tsv");
	myfile << "CloneID\tCount\n";
	for (std::pair<int,int>pr : st) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
		//myfile << "id." << pr.first << "x\t" << pr.second << "\n";
		myfile << "id." << buffer << "x\t" << pr.second << "\n";
	}
	myfile.close();
	return((int)st.size());
}
//************************************************************************************
//
//
//
//
//
//
//
//************************************************************************************
void Simulation::doSampleBlood()
{
	char buffer[100];
	sprintf(buffer,"%04i",day);
	string sDate(buffer);
	std::ofstream myfile;
	//
	// Units are used for sampling
	// Each unit represents a large but unifor group of clones
	//
	vector<int> granUnits;
	vector<int> monoUnits;
	vector<int> bUnits;
	vector<int> tUnits;
	vector<int> nkUnits;
	vector<int> wbcUnits;
	
	
	//
	// maps are used to summaraze total blood composition
	//
	std::map<int,double> mapGran;
	std::map<int,double> mapMono;
    std::map<int,double> mapWBC;
	std::map<int,double> mapB;
	std::map<int,double> mapT;
	std::map<int,double> mapNk;
	
	std::cout << "    Blood Sample: expand clones " << bloodCells.size() << "\n";
	for (auto pCell:bloodCells) {
        //
        // wbc
        //
		wbcUnits.push_back(pCell->_cloneID);
        int id = pCell->_cloneID;
        if (mapWBC.find(id) == mapWBC.end() ) {
            mapWBC[id] = (double)pCell->_count;
        } else {
            mapWBC[id] += (double)pCell->_count;
        }
        //
        // lineage
        //
		if (pCell->_type == GRAN_TYPE) {
			int id = pCell->_cloneID;
			granUnits.push_back(id);
			if (mapGran.find(id) == mapGran.end() ) {
				mapGran[id] = (double)pCell->_count;
			} else {
				mapGran[id] += (double)pCell->_count;
			}
		} else if (pCell->_type == MONO_TYPE) {
			int id = pCell->_cloneID;
			monoUnits.push_back(id);
			if (mapMono.find(id) == mapMono.end() ) {
				mapMono[id] = (double)pCell->_count;
			} else {
				mapMono[id] += (double)pCell->_count;
			}
		} else if (pCell->_type == B_TYPE) {
			int id = pCell->_cloneID;
			bUnits.push_back(id);
			if (mapB.find(id) == mapB.end() ) {
				mapB[id] = (double)pCell->_count;
			} else {
				mapB[id] += (double)pCell->_count;
			}
		} else if (pCell->_type == T_TYPE) {
			int id = pCell->_cloneID;
			tUnits.push_back(id);
			if (mapT.find(id) == mapB.end() ) {
				mapT[id] = (double)pCell->_count;
			} else {
				mapT[id] += (double)pCell->_count;
			}
		} else if (pCell->_type == NK_TYPE) {
			int id = pCell->_cloneID;
			nkUnits.push_back(id);
			if (mapNk.find(id) == mapB.end() ) {
				mapNk[id] = (double)pCell->_count;
			} else {
				mapNk[id] += (double)pCell->_count;
			}
		}
	}
	//
	// summary
	//
//	std::cout << "    Total WBC  Units         =  " << wbcUnits.size() << "\n";
//
//	std::cout << "    Total Gran Units         =  " << granUnits.size() << "\n";
//	std::cout << "    Total Gran Unique clones =  " << mapGran.size() << "\n";
//
//	std::cout << "    Total Mono Units         =  " << monoUnits.size() << "\n";
//	std::cout << "    Total Mono Unique clones =  " << mapMono.size() << "\n";
//
//	std::cout << "    Total Bcell Units         =  " << bUnits.size() << "\n";
//	std::cout << "    Total Bcell Unique clones =  " << mapB.size() << "\n";
//
//	std::cout << "    Total Tcell Units         =  " << tUnits.size() << "\n";
//	std::cout << "    Total Tcell Unique clones =  " << mapT.size() << "\n";
//
//	std::cout << "    Total NKcell Units         =  " << nkUnits.size() << "\n";
//	std::cout << "    Total NKcell Unique clones =  " << mapNk.size() << "\n";
//
//	std::cout << "       Write out acutal differentiated cell counts\n";
    //
    // wbc actual output
    //
    std::string fname = _tp->s + "Sim_" + sDate;
    myfile.open (fname + "_wbc_actual.tsv");
    myfile << "CloneID\tCount\n";
    
    for (auto pr:mapWBC) {
        char buffer[255];
        sprintf(buffer, "%06i",pr.first);
        myfile << "id." << buffer << "x" << "\t" << pr.second << endl;
    }
    myfile.close();

	
	//
	// gran actual output
	//
	fname = _tp->s + "Sim_" + sDate;
	myfile.open (fname + "_gran_actual.tsv");
	myfile << "CloneID\tCount\n";
	
	for (auto pr:mapGran) {
        char buffer[255];
        sprintf(buffer, "%06i",pr.first);
		myfile << "id." << buffer << "x" << "\t" << pr.second << endl;
	}
	myfile.close();
	//
	// mono actual output
	//

	fname = _tp->s + "Sim_" + sDate;
	myfile.open (fname + "_mono_actual.tsv");
	myfile << "CloneID\tCount\n";
	
	for (auto pr:mapMono) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
		myfile << "id." << buffer << "x" << "\t" << pr.second << endl;
	}
	myfile.close();
	//
	// b actual output
	//
	fname = _tp->s + "Sim_" + sDate;
	myfile.open (fname + "_b_actual.tsv");
	myfile << "CloneID\tCount\n";
	
	for (auto pr:mapB) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
		myfile << "id." << buffer << "x" << "\t" << pr.second << endl;
	}
	myfile.close();
	//
	// t actual output
	//
	fname = _tp->s + "Sim_" + sDate;
	myfile.open (fname + "_t_actual.tsv");
	myfile << "CloneID\tCount\n";
	
	for (auto pr:mapT) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
		myfile << "id." << buffer << "x" << "\t" << pr.second << endl;
	}
	myfile.close();
	//
	// nk actual output
	//
	fname = _tp->s + "Sim_" + sDate;
	myfile.open (fname + "_nk_actual.tsv");
	myfile << "CloneID\tCount\n";
	
	for (auto pr:mapNk) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
		myfile << "id." << buffer << "x" << "\t" << pr.second << endl;
	}
	myfile.close();
	//
	// sampling
	//
	std::cout << "       Sample Lineages\n";
	doSampleAll(mapWBC,"WBC_GFP+",sDate,_tp->s);
	doSampleAll(mapGran,"Gran",sDate,_tp->s);
	
	doSample(wbcUnits,"i_WBC_GFP+",sDate,_tp->s);
	doSample(granUnits,"i_Gran",sDate,_tp->s);
	doSample(monoUnits,"Mono",sDate,_tp->s);
	doSample(bUnits,"Bcell",sDate,_tp->s);
	doSample(tUnits,"Tcell",sDate,_tp->s);
	doSample(nkUnits,"NKcell",sDate,_tp->s);
	
}
//************************************************************************************
//
// Sim Blood Cells
//
//
//  MPP immediately go into one lineage
//
//
//************************************************************************************
void Simulation::simBloodCells()
{
	std::vector<BloodCell*> bloodCellsTemp;
	std::vector<BloodCell*> deadCells;
	long bloodGranCount = 0;
	long bloodMonoCount = 0;
	long bloodBCellCount = 0;
	long bloodTCellCount = 0;
	long bloodNKCellCount = 0;
	for (int i = 0; i < bloodCells.size(); i++) {
		BloodCell *pCell = bloodCells[i];
	
		if (pCell->_cloneID > _tp->cells) {
			cout << "Bad clone ID " << pCell->_cloneID << endl;
		}
		//
		// need stochastic cell removal
		//
		if (pCell->_type == GRAN_TYPE) {
			pCell->_age++;
			if (pCell->_age <= 1) {
				bloodGranCount += pCell->_count;
				bloodCellsTemp.push_back(pCell);
			} else {
				pCell->_cloneID = 1000000;
				deadCells.push_back(pCell);
			}
		} else if (pCell->_type == MONO_TYPE) {
			pCell->_age++;
			// 7
			if (pCell->_age <= 2) {
				bloodMonoCount += pCell->_count;
				bloodCellsTemp.push_back(pCell);
			} else {
				pCell->_cloneID = 1000000;
				deadCells.push_back(pCell);
			}
		} else if (pCell->_type == B_TYPE) {
			pCell->_age++;
			//
			if (pCell->_age <= 28) {
				bloodBCellCount += pCell->_count;
				bloodCellsTemp.push_back(pCell);
			} else {
				pCell->_cloneID = 1000000;
				deadCells.push_back(pCell);
			}
		} else if (pCell->_type == T_TYPE) {
			pCell->_age++;
			//
			if (pCell->_age <= 56) {
				bloodTCellCount += pCell->_count;
				bloodCellsTemp.push_back(pCell);
			} else {
				pCell->_cloneID = 1000000;
				deadCells.push_back(pCell);
			}
		} else if (pCell->_type == NK_TYPE) {
			pCell->_age++;
			//
			if (pCell->_age <= 14) {
				bloodNKCellCount += pCell->_count;
				bloodCellsTemp.push_back(pCell);
			} else {
				pCell->_cloneID = 1000000;
				deadCells.push_back(pCell);
			}
		} else {
			std::cout << "Error bad cell type " << pCell->_type << "\n";
		}
	}
	
	
//	std::cout << "    blood cell clones = " << bloodCellsTemp.size() << endl;
//	std::cout << "                blood gran cell count = " << bloodGranCount << "\n";
//	std::cout << "                blood mono cell count = " << bloodMonoCount << "\n";
//	std::cout << "                blood B    cell count = " << bloodBCellCount << "\n";
//	std::cout << "                blood T    cell count = " << bloodTCellCount << "\n";
//	std::cout << "                blood NK   cell count = " << bloodNKCellCount << "\n";
	
	bloodCells.clear();
	for (auto pCell:bloodCellsTemp) {
		bloodCells.push_back(pCell);
	}
	
	for (auto p:bloodCells) {
		if (p->_cloneID > _tp->cells) {
			cout << "Bad clone ID " << p->_cloneID << endl;
		}
	}
	bloodCellsTemp.clear();
	for (int i = 0; i < deadCells.size();i++) {
		delete deadCells[i];
	}
}



//************************************************************************************
//
// Sim MPP Cells
//
//
//  MPP immediately go into one lineage
//
//
//************************************************************************************
void Simulation::simMppCellsV1()
{
	int newGran  = 0;
	int newMono  = 0;
	int newB     = 0;
	int newT     = 0;
	int newNK    = 0;
	int ery      = 0;
	int mppDie   = 0;
	int mppExpansion = 1;
	
	std::cout << "    Starting Total MPP cells before division = " << mppCells.size() << "\n";
	
	
	for (auto pCell:mppCells) {
		for (int i=0;i<mppExpansion;i++) {
			//
			// mpp die?
			//
			int rmpp = rand() % 100;
			if (rmpp < 98) {
				int r2 = rand() % 100;
				if (r2 <= 33) {
					// erythroid: nothing happens
					ery++;
				} else if (r2 <= 66) {
					int r3 = rand() % 100;
					if (r3 <= 50) {
						GranCell *pGran = new GranCell(pCell->_cloneID,0,pCell->_age,pCell->_ef);
						granCells.push_back(pGran);
						newGran++;
					} else {
						MonoCell *pMono = new MonoCell(pCell->_cloneID,0,pCell->_age,pCell->_ef);
						monoCells.push_back(pMono);
						newMono++;
					}
				} else {
					int r3 = rand() % 100;
					if (r3 <= 33) {
						BCell *pB = new BCell(pCell->_cloneID,0,pCell->_age,pCell->_ef);
						bCells.push_back(pB);
						newB++;
					} else if (r3 <= 66) {
						TCell *pT = new TCell(pCell->_cloneID,0,pCell->_age,pCell->_ef);
						tCells.push_back(pT);
						newT++;
					} else {
						NkCell *pNK = new NkCell(pCell->_cloneID,0,pCell->_age,pCell->_ef);
						nkCells.push_back(pNK);
						newNK++;
					}
				}
			} else {
				mppDie++;
			}
		}
		// mpp is gone
		pCell->_cloneID = 1000000;
		delete pCell;
	}
	mppCells.clear();
	std::cout << "      Ery from MPP      = " << ery << "\n";
	std::cout << "      New Gran From MPP = " << newGran << "\n";
	std::cout << "      New Mono From MPP = " << newMono << "\n";
	std::cout << "      New B From MPP    = " << newB << "\n";
	std::cout << "      New T From MPP    = " << newT << "\n";
	std::cout << "      New NK From MPP   = " << newNK << "\n";
	std::cout << "      MPP die           = " << mppDie << "\n";
}


//************************************************************************************
//
// V1: assymetric differentiation
//
//
//
//
//
//************************************************************************************
void Simulation::simStemCellsV1()
{
	//
	// dividing stem cells progress along cycle
	//
	std::cout << "    cycling stem cells = " << stemCellsCycling.size() << "\n";
	std::vector<StemCell*> tContinueCycle;
	//std::vector<StemCell*> tNewResting;
	std::vector<StemCell*> tNewRecovering;
	std::vector<StemCell*> tSurvieApop;
	std::vector<StemCell*> result;
	std::vector<StemCell*> deadCells;
	//
	// stem cells finish recovery
	//
	for (int i = 0; i < stemCellsRecovering.size(); i++) {
		StemCell* pCell = stemCellsRecovering[i];
		pCell->_recoveryCount++;
		if (pCell->_recoveryCount >= pCell->_recoveryTime) {
			stemCellsResting.push_back(pCell);
		} else {
			tNewRecovering.push_back(pCell);
		}
	}
	stemCellsRecovering.clear();
	for (auto pCell:tNewRecovering) {
		stemCellsRecovering.push_back(pCell);
	}
	tNewRecovering.clear();
	//
	// cells that finish cycle
	//
	int newCells = 0;
	int diffCount = 0;
	int stemCount = 0;
	for (int i = 0; i < stemCellsCycling.size(); i++) {
		StemCell* pCell = stemCellsCycling[i];
		pCell->_cycleCount ++;
		if (pCell->_cycleCount >= 1) {
			//cout << pCell->_recoveryTime << "\t    " <<  pCell->_epigeneticFactor << "   " << pCell->_cloneID << "\n";
			newCells += 2;
			//
			// cell cycle completes: sym or asym division
			//   chance stem = 0.5
			//   chance diff = 0.5
			//
			pCell->_cycleCount = 0;
			pCell->_age++;
			//
			// for each daughter...stem or mpp
			//  for stem is there a recovery time or not?
			//
			long it = rand();
			long id = RAND_MAX;
			double r2 = (double)it / (double)id;
			if (r2 < pCell->_pStem) {
				
				if (pCell->_cloneID > _tp->cells) {
					std::cout << "ERROR bad cloneID\n";
				}
				
				StemCell *pNewDaughter = new StemCell(pCell->_cloneID,pCell->_age,pCell->_epigeneticFactor,pCell->_recoveryTime);
				pNewDaughter->_pStem = pCell->_pStem;
				pNewDaughter->_pDif = pCell->_pDif;
				pNewDaughter->_pApop = pCell->_pApop;
				pNewDaughter->_pCycle = pCell->_pCycle;
				
				if (pNewDaughter->_recoveryTime == 0) {
					stemCellsResting.push_back(pNewDaughter);
				} else {
					stemCellsRecovering.push_back(pNewDaughter);
				}
				stemCount++;
			} else if (r2 < (pCell->_pStem + pCell->_pDif)) {
				MppCell *pMpp = new MppCell(pCell->_cloneID,pCell->_age,pCell->_pCycle);
				mppCells.push_back(pMpp);
				diffCount++;
			} else {
				// cell dies
			}
			//
			// do this one last so pCell not deleted too soon!
			//
			it = rand();
			id = RAND_MAX;
			double r1 = (double)it / (double)id;
			
			if (r1 < pCell->_pStem) {
				pCell->_recoveryCount = 0;
				if (pCell->_recoveryTime == 0) {
					stemCellsResting.push_back(pCell);
				} else {
					stemCellsRecovering.push_back(pCell);
				}
				stemCount++;
			} else if (r1 < (pCell->_pStem + pCell->_pDif)){
				MppCell *pMpp = new MppCell(pCell->_cloneID,pCell->_age,pCell->_pCycle);
				mppCells.push_back(pMpp);
				pCell->_cloneID = 1000000;
				pCell->_state = 411;
				delete pCell;
				diffCount++;
			} else {
				// cell dies
			}
		} else {
			tContinueCycle.push_back(pCell);
		}
	}
	std::cout << "    stem cells completing cell cycle = " << newCells << "\n";
	std::cout << "      New Stem count = " << stemCount << "\n";
	std::cout << "      New Diff count = " << diffCount << "\n";
	
	stemCellsCycling.clear();
	for (auto pCell:tContinueCycle) {
		stemCellsCycling.push_back(pCell);
	}
	tContinueCycle.clear();
//	//
//	// add cells completing cycle back to resting pool
//	//
//	for (auto pCell:tNewResting) {
//		stemCellsResting.push_back(pCell);
//	}
	//
	// some stem cells die
	//
	//if (stemCellsResting.size() > 0) {
	
	if ((long)stemCellsResting.size() > (long)_apopLimit) {
		std::cout << "Over pop limit   " << stemCellsResting.size() << "\n";
		//int nDie = (int)(stemCellsResting.size() * _tp->deathRate);
		int nDie = (int)(stemCellsResting.size() - _apopLimit);
		std::sample(stemCellsResting.begin(),
					stemCellsResting.end(),
					std::back_inserter(deadCells),
					nDie,
					std::mt19937{std::random_device{}()});
		
		cout << "    cells die from apop = " << deadCells.size() << "\n";
		for (auto pCell:deadCells) {
			pCell->_alive = 0;
		}
	}
	//
	// remake,except for dead cells
	//
	for (auto pCell:stemCellsResting) {
		if (pCell->_alive == 1) {
			tSurvieApop.push_back(pCell);
		}
	}
	stemCellsResting.clear();
	//
	// a fraction of resting stem cells start to divide
	//
	result.clear();
	//double repRate = _tp->r1;
	//int rep = (int)floor(tSurvieApop.size() * repRate);
	//
	// sample from vector
	//
	//
	//	std::sample(tSurvieApop.begin(),
	//				tSurvieApop.end(),
	//				std::back_inserter(result),
	//				rep,
	//				std::mt19937{std::random_device{}()});
	//
	// select cycle cells
	//
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	int newCycle = 0;
	for (auto pCell:tSurvieApop) {
		double r = pCell->_pCycle * distribution(generator);
		if (r >= _tp->r1) {
			// cycleCount = 1 will cause cell not to be added back to stemCellResting
			pCell->_cycleCount = 1;
			stemCellsCycling.push_back(pCell);
			newCycle++;
		}
	}

	//
	// move survive apop cells to resting list,
	// except for newly cycling cells
	//
	for (auto pCell:tSurvieApop) {
		if (pCell->_cycleCount == 0) {
			stemCellsResting.push_back(pCell);
		}
	}
	tSurvieApop.clear();
	//
	// move cells completing cell cycle to list
	//
	long tc = stemCellsResting.size() + stemCellsCycling.size() + stemCellsRecovering.size();
	std::cout << "    New cycle cells             = " << newCycle << "\n";
	std::cout << "    resting cells               = " << stemCellsResting.size() << "\n";
	std::cout << "    total stem cells recovering = " << stemCellsRecovering.size() << "\n";
	std::cout << "    total stem cells            = " << tc << "\n";

//	if (day == 10) {
//		for (auto pCell:stemCellsRecovering) {
//			cout << pCell->_recoveryCount << "\t" << pCell->_recoveryTime << "\n";
//		}
//	}
	
	for (auto p:stemCellsCycling) {
		if (p->_cloneID > _tp->cells*10) {

			std::cout << "ERROR bad id \n";
		}
	}
	for (auto p:stemCellsResting) {
		if (p->_cloneID > _tp->cells*10) {
			std::cout << "ERROR bad id \n";
		}
	}
//	for(auto p:deadCells) {
//		delete p;
//	}
//	deadCells.clear();
}



//************************************************************************************
//
// V1: assymetric differentiation
//
//
//
//
//
//************************************************************************************
Simulation::Simulation(TEST_PARAM* ptp)
{
	_tp = ptp;
	bfGran = new BloodFIFO();
	bfMono = new BloodFIFO();
	bfB = new BloodFIFO();
	bfT = new BloodFIFO();
	bfNK = new BloodFIFO();
}
//************************************************************************************
//
//
//
//
//
//
//
//************************************************************************************
void Simulation::run_sim() {
	std::string outputDirectory = _tp->s;

	
	int sampleDays[] = {0,10,20,30,40,50,60,70,80,90,
		100,150,200,250,300,350,400,450,500,
		600,700,800,900,1000,1100,1200,1300,1400,1500};
	//
	// init stem cells
	//
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(_tp->epiMean,_tp->epiVar);
	
	for (int i = 0; i <_tp->cells;i++) {
		StemCell *pst = new StemCell(i,0,0,0);
		//
		//	stem cells have a epigenetic cycle
		//  time of a normal distribution around X days
		//  divided into 6
		//
		double number = distribution(generator);
		if (number < 0.1) number = 0.1;
		//if (number > 1.0) number = 1.0;
		//std::cout << "             " << number << "\n";
		pst->_pCycle = number;
		pst->_pStem = _tp->pStem;
		pst->_pDif = 1.0-pst->_pStem;
		pst->_pApop = 0.0;
//		int population = rand() % 100;
//
//		if (population < 50) {
//			pst->_pStem = 0.51;
//			pst->_pDif  = 0.49;
//			pst->_pApop = 0.0;
//		} else {
//			pst->_pStem = 0.40;
//			pst->_pDif  = 0.6;
//			pst->_pApop = 0.0;
//		}
		pst->_epigeneticFactor = 1.0;
		pst->_recoveryTime = (int)0;
		stemCellsResting.push_back(pst);
	}
	
	//
	// stem cell details
	//
	std::ofstream myfile;
	myfile.open (_tp->s + "initial_stem_id.tsv");
	myfile << "CloneID\tCount\n";
	for (auto pCell:stemCellsResting) {
		char buffer[255];
		sprintf(buffer, "%06i",pCell->_cloneID);
		myfile << "id." << buffer << "x\t" << pCell->_epigeneticFactor << "\t  " << pCell->_recoveryTime << "\n";
	}
	myfile.close();
    //
    // write out globalID
    //
    // ID    count    NormCount    refseq    multi    align
    // Z15086_000001    12256    2107767    GATTTAACAAAAGTTCATGAAAAATAACTT    S    chr7_18150699_-
    std::ofstream fgl;
    std::string filename = outputDirectory + "globalID.txt";
    fgl.open (filename);
    fgl << "ID\tcount\tNormCount\trefseq\tmulti\talign\n";
    
    for (int i = 0; i < _tp->cells;i++) {
        char buffer[255];
        //sprintf(buffer, "ZSIM_%06i",i);
        sprintf(buffer, "id.%06ix",i);
        
        fgl << buffer << "\t" << 1 << "\t" << 1 << "\t" << "_\t" << "_\t" << "_\n";
    }
    fgl.close();
	//
	// sim
	//
	for (day = 0; day <= 1500; day++) {
		cout << "DAY  " << day << " -----------------------------------" << endl;
		//
		// simulate blood cells
		//
		simBloodCells();
		//
		// get newly matured gran cells and add to blood
		//
		bfGran->get(&bloodCells);
		//
		// get newly matured mono cells
		//
		bfMono->get(&bloodCells);
		//
		// get newly matured b cells
		//
		bfB->get(&bloodCells);
		//
		// get newly matured t cells
		//
		bfT->get(&bloodCells);
		//
		// get newly matured nk cells
		//
		bfNK->get(&bloodCells);
		//
		// sim mpp cells
		//
        simMppCellsV1();
		//
		//
		//
		cout << "    add " << granCells.size() << " Gran clones to schedule fifo\n";
		//cout << "    add " << monoCells.size() << " Mono clones to schedule fifo\n";
		//cout << "    add " << bCells.size() << " B Cell clones to schedule fifo\n";
		cout << "    add " << tCells.size() << " T Cell clones to schedule fifo\n";
		//cout << "    add " << nkCells.size() << " NK Cell clones to schedule fifo\n";
		//
		//  add new grans to development schedule
		//
		// expansion factor = days expanding (makes 2 ^ expansionFactor cells)
		//
		int expansionFactor = 10;
		int matureDelay = 14;
		bfGran->distribute(&granCells,matureDelay,expansionFactor,2.0, 0.0);
		//bfGran->Diag();
		//
		//  add new monos to development schedule
		//  
		if (day < 50) {
			monoCells.clear();
			bCells.clear();
			tCells.clear();
			nkCells.clear();
		} else {
			expansionFactor = 4;
			matureDelay = 16;
			bfMono->distribute(&monoCells,matureDelay,expansionFactor, 2.0, 0.0);
			//bfMono->Diag();
			//
			// add new bCells to development schedule
			//
			expansionFactor = 4;
			matureDelay = 20;
			bfB->distribute(&bCells,matureDelay,expansionFactor, 2.0, 2.0);
			//
			// add new bCells to development schedule
			//
			expansionFactor = 4;
			matureDelay = 40;
			bfT->distribute(&tCells,matureDelay,expansionFactor, 2.0, 2.0);
			//
			// add new nkCells to['09']development schedule
			//
			expansionFactor = 4;
			matureDelay = 20;
			bfNK->distribute(&nkCells,matureDelay,expansionFactor, 2.0, 2.0);
		}


		//
		// simulate cycling and differentiation of stem cells
		//
		
		simStemCellsV1();
		
//		cout << "    add " << mppCells.size()-old << " new mpp cells, total = " << mppCells.size() << "\n";

		//***********************************************************
		//
		// END of simulation
		//
		// keep track of unique clone count
		//***********************************************************
		bool match = false;
		for (int sd:sampleDays) {
			if (day == sd) {
				match = true;
				break;
			}
		}
		if (match) {
			int unique = doSampleStem();
            std::pair<int,int> pr(day,unique);
			uniqueCloneCount.push_back(pr);
			doSampleBlood();
		}
		//
		// switch date params if needed   !!! yes this happens each loop
		//
		if (day > _tp->day) {
			_tp->r1 = _tp->r2;
		}
	}

	//
	// write out uniques
	//
	std::ofstream mySfile;
	filename = outputDirectory + "uniqueCloneCount.tsv";
	mySfile.open (filename);
	mySfile << "day\tCount\n";
	int day = 0;
    for (auto pr:uniqueCloneCount) {
		mySfile << pr.first << "\t" << pr.second << "\n";
		day += 10;
	}
	mySfile.close();
}
