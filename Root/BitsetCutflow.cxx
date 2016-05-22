#include <MetTriggerPackage/BitsetCutflow.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(BitsetCutflow)

BitsetCutflow::BitsetCutflow(EL::Worker *wk){
  m_wk = wk;
  m_cutflowHist = new TH1I("cutflow_hist","Cutflow",99,-0.5,98.5);
  m_wk->addOutput(m_cutflowHist);
  m_cutflowStepCounter = 1;
  m_cutflowHist->GetXaxis()->SetBinLabel(m_cutflowStepCounter,"All");
  m_cutflowHist->Fill(0.0,-1.0); /// to compensate for the first event false filling
}

BitsetCutflow::~BitsetCutflow(){
  
}

void BitsetCutflow::PushBitSet(){
  
  /// First bin: All events
  m_cutflowHist->Fill(0);
  
  /// first and second bits will be always 0
  for (unsigned i=1; i<=m_mapCutflowSteps.size(); i++){
    /// first bit is always zero
    if (m_bitset[i+1]!=m_bitset[0]) m_cutflowHist->Fill(i); 
  }
  
  /// clean bit for next event
  m_bitset.reset();
}

void BitsetCutflow::FillCutflow(string cutflowStepName, bool fillCutflow){
  /// check if there is a slot already created for this cutflow step
  if (fillCutflow){
    unsigned int slotPosition = m_mapCutflowSteps[cutflowStepName];
    if (slotPosition==0){
      m_cutflowStepCounter++;
      m_mapCutflowSteps[cutflowStepName] = m_cutflowStepCounter;
      slotPosition = m_cutflowStepCounter;
      m_cutflowHist->GetXaxis()->SetBinLabel(m_cutflowStepCounter,cutflowStepName.c_str());
    }
    m_bitset[slotPosition]=1;
  }
}

void BitsetCutflow::PrintCutflowLocally(){
  for (int i=1; i<=m_cutflowHist->GetNbinsX(); i++){
    int binContent = m_cutflowHist->GetBinContent(i);
    if (binContent<=0) break;
    string binLabel = m_cutflowHist->GetXaxis()->GetBinLabel(i);
    cout << binLabel << ":\t" << binContent << endl;
  }
}
