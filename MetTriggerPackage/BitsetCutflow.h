#ifndef BitsetCutflow_H
#define BitsetCutflow_H

#include <TH1I.h>
#include <iostream>
#include <map>
#include <bitset>

#include "EventLoop/Worker.h"

using namespace std;

class BitsetCutflow
{

public:
	BitsetCutflow(EL::Worker *wk);
	~BitsetCutflow();

	void FillCutflow(string cutflowStepName, bool fillCutflow = true);
	void PrintCutflowLocally();
        
	/// WARNING call this function on the BEGIN of EVENT!!!
	/// WARNING call this function in the finalize() function!!!
	void PushBitSet();
	
private:
	
	unsigned int m_cutflowStepCounter; //!
	
	/// link to EventLoop worker;
	EL::Worker *m_wk; //!
	
	/// actual bitset were we will store cutflow on event-by-event basis
	bitset<99> m_bitset; //!
	
	/// map with cutflow steps
	map<string,unsigned int> m_mapCutflowSteps; //!
	
	TH1I* m_cutflowHist; //!
	
	/// this is needed to distribute the algorithm to the workers
	ClassDef(BitsetCutflow, 1);

};

#endif
