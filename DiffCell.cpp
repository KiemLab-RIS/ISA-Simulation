//
//  DiffCell.cpp
//  StemSim2
//
//  Created by mark enstrom on 4/16/21.
//  Copyright © 2021 Mark Enstrom. All rights reserved.
//
#include "unordered_set"
#include "DiffCell.hpp"

std::unordered_set<long> BobFloydAlgo(long sampleSize, long rangeUpperBound)
{
     std::unordered_set<long> sample;
     std::default_random_engine generator;

     for(long d = rangeUpperBound - sampleSize; d < rangeUpperBound; d++)
     {
		 double f = std::uniform_real_distribution<>(0, 1.0)(generator);
		 long t = (long)floor(f*rangeUpperBound);
		 //long t = uniform_int_distribution<>(0, d)(generator);
		 if (sample.find(t) == sample.end() )
		   sample.insert(t);
		 else
		   sample.insert(d);
     }
     return sample;
}

struct myclassDifSort {
    bool operator() (std::pair<int,int> pa1,std::pair<int,int> pa2) { return (pa1.second > pa2.second);}
} DifSortCompObj;


void doSampleAll(std::map<int,double> &cloneCounts,
						std::string cellType,
						std::string sDate,
                        std::string outputDirectory)
{
	long totalCloneCount = 0;
	std::vector<std::pair<int,int>> sam;
	
	
	for (auto pr:cloneCounts) {
		totalCloneCount += pr.second;
	}
	
	if (totalCloneCount > 0) {
		
		std::cout << "total clone count for sampling = " << totalCloneCount << "\n";

		std::vector<int> idVec(totalCloneCount,0);
		
		int index = 0;
		for (auto pr:cloneCounts) {
			for (int j =0; j < pr.second;j++) {
				idVec[index++] = pr.first;
			}
		}
		
		//
		//long sampleSize = 25000;// (500,000 cells * 0.05)
		
		long sampleSize = 20000;// (500,000 cells * 0.05)
		
		
		std::unordered_set st = BobFloydAlgo(sampleSize, totalCloneCount);
		
		std::map<int,int> mapSampleG;
		
		for (auto l:st) {
			int cloneID = idVec[l];
			
			if (mapSampleG.find(cloneID) == mapSampleG.end() ) {
				mapSampleG[cloneID] = 1;
			} else {
				mapSampleG[cloneID] ++;
			}
		}
		//
		// make into vector
		//

		
		for (auto pr:mapSampleG) {
			sam.push_back(pr);
		}
		//
		// sort
		//
		std::sort(sam.begin(),sam.end(),DifSortCompObj);
		
	}
	//
	//
	//
	std::ofstream myfile;

    //   Z15086_0070DPT_PB_Gran.tsv

	std::string fname = outputDirectory + "ZSIM_" + sDate + "DPT";
	fname = fname + "_PB_" + cellType + ".tsv";

	myfile.open (fname);

    myfile << "id\tGlobalIndex\toCount\tnCount\tfrac\tlog2\tfracLog2\tspan\tmultiR\talign\trefseq\n";

    double totalCount = 0.0;
    for (auto pr:sam) {
        totalCount += pr.second;
    }

	for (auto pr:sam) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
        myfile << "id." << buffer << "x" << "\t";
        myfile << pr.first << "\t";
        myfile << pr.second << "\t";
        double nCount = pr.second/totalCount * 1000000.0;
        myfile << (int)floor(nCount) << "\t";
        myfile << nCount/1000000.0 << "\t";
        myfile << 0.0 << "\t";
        myfile << 0.0 << "\t";
        myfile << 0 << "_\t";
        myfile << "_\t";
        myfile << "_\t";
        myfile << "_" << std::endl;
	}
	myfile.close();

}





void doSample(std::vector<int> cellUnits,
						std::string cellType,
						std::string sDate,
                        std::string outputDirectory)
{
	if (cellUnits.size() <= 0) return;

	int sampleSize = (int)500000.0 * 0.05;
	//
	// collect random samples
	//
	std::map<int,int> mapSampleG;
	std::vector<int> sam;
	std::sample(cellUnits.begin(),
				cellUnits.end(),
				std::back_inserter(sam),
				sampleSize,
				std::mt19937{std::random_device{}()});
	//
	// small pcr sim
	//
	std::vector<int> pcr;
	for (int i = 0; i < 4; i++) {
		for (int id:sam) {
			int r2 = rand() % 100;
			if (r2 <= 80) {
				pcr.push_back(id);
			}
		}
		for (int id:pcr) {
			sam.push_back(id);
		}
		pcr.clear();
	}
	//
	// simulate sequencing
	//
	int risCount = 10000;
	std::vector<int> ris;
	std::sample(sam.begin(),
				sam.end(),
				std::back_inserter(ris),
				risCount,
				std::mt19937{std::random_device{}()});
	//
	// build a map to summarize sequnceing results
	//
	for (int sampleID:ris) {
		if (mapSampleG.find(sampleID) == mapSampleG.end() ) {
			mapSampleG[sampleID] = 1;
		} else {
			mapSampleG[sampleID] += 1;
		}
	}
	//
	//
	//
	std::ofstream myfile;
	
    //   Z15086_0070DPT_PB_Gran.tsv
    
	std::string fname = outputDirectory + "ZSIM_" + sDate + "DPT";
	fname = fname + "_PB_" + cellType + ".tsv";
	
	myfile.open (fname);
	
    myfile << "id\tGlobalIndex\toCount\tnCount\tfrac\tlog2\tfracLog2\tspan\tmultiR\talign\trefseq\n";
    
    double totalCount = 0.0;
    for (auto pr:mapSampleG) {
        totalCount += pr.second;
    }
	
	for (auto pr:mapSampleG) {
		char buffer[255];
		sprintf(buffer, "%06i",pr.first);
        myfile << "id." << buffer << "x" << "\t";
        myfile << pr.first << "\t";
        myfile << pr.second << "\t";
        double nCount = pr.second/totalCount * 1000000.0;
        myfile << (int)floor(nCount) << "\t";
        myfile << nCount/1000000.0 << "\t";
        myfile << 0.0 << "\t";
        myfile << 0.0 << "\t";
        myfile << 0 << "_\t";
        myfile << "_\t";
        myfile << "_\t";
        myfile << "_" << std::endl;
	}
	myfile.close();

}
