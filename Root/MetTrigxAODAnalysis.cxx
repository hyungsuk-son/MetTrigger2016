#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MetTriggerPackage/MetTrigxAODAnalysis.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODBase/IParticleHelpers.h"
#include "AthContainers/ConstDataVector.h"

#include <TSystem.h>
#include <TFile.h>

#include "xAODRootAccess/tools/Message.h"

#include "PATInterfaces/CorrectionCode.h" // to check the return correction code status of tools
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/AuxContainerBase.h"

static std::string jetType = "AntiKt4EMTopoJets";

// Global accessors and decorators
static SG::AuxElement::Decorator<char> dec_baseline("baseline");
static SG::AuxElement::Decorator<char> dec_signal("signal");
static SG::AuxElement::Decorator<char> dec_bad("bad");
static SG::AuxElement::Accessor<float>  acc_jvt("Jvt");
static SG::AuxElement::ConstAccessor<float> cacc_jvt("Jvt");
// For ORTools
static const std::string inputLabel = "selected";
static const std::string outputLabel = "overlaps";
const bool outputPassValue = false;
static const std::string bJetLabel = "";
//static SG::AuxElement::Accessor<char> overlapAcc("overlaps");
ort::inputAccessor_t selectAcc(inputLabel);
ort::inputDecorator_t selectDec(inputLabel);
ort::outputAccessor_t overlapAcc(outputLabel);
ort::inputDecorator_t bJetDec(bJetLabel);
ort::objLinkAccessor_t objLinkAcc("overlapObject");
//static const bool outputPassValue = false;
//static const std::string outputLabel = outputPassValue? "passOR" : "overlaps";
//static SG::AuxElement::Decorator<char> dec_overlap(outputLabel);

// Scale Factor decorators
static SG::AuxElement::Decorator<double> dec_scalefactor("scalefactor");

struct DescendingPt:std::function<bool(const xAOD::IParticle*, const xAOD::IParticle*)> {
  bool operator()(const xAOD::IParticle* l, const xAOD::IParticle* r)  const {
    return l->pt() > r->pt();
  }
};


// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
  do {                                                     \
    if( ! EXP.isSuccess() ) {                             \
      Error( CONTEXT,                                    \
          XAOD_MESSAGE( "Failed to execute: %s" ),    \
#EXP );                                     \
      return EL::StatusCode::FAILURE;                    \
    }                                                     \
  } while( false )



// this is needed to distribute the algorithm to the workers
ClassImp(MetTrigxAODAnalysis)



MetTrigxAODAnalysis :: MetTrigxAODAnalysis ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode MetTrigxAODAnalysis :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  job.useXAOD ();
  xAOD::Init().ignore(); // call before opening first file
  //  CP::CorrectionCode::enableFailure();
  EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.


  // Average and Actual Interaction
  h_avg_interaction = new TH1F("h_avg_interaction", "Average Interactions per Crossing;<#mu>", 40, 0, 40);
  h_act_interaction = new TH1F("h_act_interaction", "Actual Interactions per Crossing;<#mu>", 40, 0, 40);
  wk()->addOutput (h_avg_interaction);
  wk()->addOutput (h_act_interaction);

  // Bunch Crossing ID (BCID)
  h_bcid = new TH1F("h_bcid", "Bunch Crossing ID;BCID", 3600, 0, 3600);
  wk()->addOutput (h_bcid);

  // L1 MET
  h_l1_mex = new TH1F("h_l1_mex", "L1 METx (GeV);METx (GeV)", 150, -150,  150); // L1 METx [GeV]
  h_l1_mey = new TH1F("h_l1_mey", "L1 METy (GeV);METy (GeV)", 150, -150,  150); // L1 METy [GeV]
  h_l1_met = new TH1F("h_l1_met", "L1 |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // L1 MET [GeV]
  h_l1_sumet = new TH1F("h_l1_sumet", "L1 SumEt (GeV);SumEt (GeV)", 250, 0, 2000); // L1 SumET [GeV]
  h_l1_phi = new TH1F("h_l1_phi", "L1 MET #phi (rad);MET #phi (rad)", 32, -3.1416, 3.1416); // L1 phi [GeV]
  wk()->addOutput (h_l1_mex);
  wk()->addOutput (h_l1_mey);
  wk()->addOutput (h_l1_met);
  wk()->addOutput (h_l1_sumet);
  wk()->addOutput (h_l1_phi);

  // HLT Cell MET
  h_hlt_ex = new TH1F("h_hlt_ex", "HLT (CELL) Missing E_{x};E_{x} (GeV)", 150, -150,  150); // HLT MEx [GeV]
  h_hlt_ey = new TH1F("h_hlt_ey", "HLT (CELL) Missing E_{y};E_{y} (GeV)", 150, -150,  150); // HLT MEy [GeV]
  h_hlt_met = new TH1F("h_hlt_met", "HLT (CELL) |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // HLT MET [GeV]
  h_hlt_sumet = new TH1F("h_hlt_sumet", "HLT (CELL) Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // HLT SumET [GeV]
  h_hlt_phi = new TH1F("h_hlt_phi", "HLT (CELL) MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // HLT phi [GeV]
  wk()->addOutput (h_hlt_ex);
  wk()->addOutput (h_hlt_ey);
  wk()->addOutput (h_hlt_met);
  wk()->addOutput (h_hlt_sumet);
  wk()->addOutput (h_hlt_phi);

  // HLT MHT MET
  h_hlt_mht_ex = new TH1F("h_hlt_mht_ex", "HLT (MHT) Missing E_{x};E_{x} (GeV)", 150, -150,  150); // HLT MEx [GeV]
  h_hlt_mht_ey = new TH1F("h_hlt_mht_ey", "HLT (MHT) Missing E_{y};E_{y} (GeV)", 150, -150,  150); // HLT MEy [GeV]
  h_hlt_mht_met = new TH1F("h_hlt_mht_met", "HLT (MHT) |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // HLT MET [GeV]
  h_hlt_mht_sumet = new TH1F("h_hlt_mht_sumet", "HLT (MHT) Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // HLT SumET [GeV]
  h_hlt_mht_phi = new TH1F("h_hlt_mht_phi", "HLT (MHT) MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // HLT phi [GeV]
  wk()->addOutput (h_hlt_mht_ex);
  wk()->addOutput (h_hlt_mht_ey);
  wk()->addOutput (h_hlt_mht_met);
  wk()->addOutput (h_hlt_mht_sumet);
  wk()->addOutput (h_hlt_mht_phi);

  // HLT Topocl MET
  h_hlt_topocl_ex = new TH1F("h_hlt_topocl_ex", "HLT (topocl) Missing E_{x};E_{x} (GeV)", 150, -150,  150); // HLT MEx [GeV]
  h_hlt_topocl_ey = new TH1F("h_hlt_topocl_ey", "HLT (topocl) Missing E_{y};E_{y} (GeV)", 150, -150,  150); // HLT MEy [GeV]
  h_hlt_topocl_met = new TH1F("h_hlt_topocl_met", "HLT (topocl) |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // HLT MET [GeV]
  h_hlt_topocl_sumet = new TH1F("h_hlt_topocl_sumet", "HLT (topocl) Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // HLT SumET [GeV]
  h_hlt_topocl_phi = new TH1F("h_hlt_topocl_phi", "HLT (topocl) MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // HLT phi [GeV]
  wk()->addOutput (h_hlt_topocl_ex);
  wk()->addOutput (h_hlt_topocl_ey);
  wk()->addOutput (h_hlt_topocl_met);
  wk()->addOutput (h_hlt_topocl_sumet);
  wk()->addOutput (h_hlt_topocl_phi);

  // HLT Topocl_ps MET
  h_hlt_topocl_ps_ex = new TH1F("h_hlt_topocl_ps_ex", "HLT (topocl_PS) Missing E_{x};E_{x} (GeV)", 150, -150,  150); // HLT MEx [GeV]
  h_hlt_topocl_ps_ey = new TH1F("h_hlt_topocl_ps_ey", "HLT (topocl_PS) Missing E_{y};E_{y} (GeV)", 150, -150,  150); // HLT MEy [GeV]
  h_hlt_topocl_ps_met = new TH1F("h_hlt_topocl_ps_met", "HLT (topocl_PS) |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // HLT MET [GeV]
  h_hlt_topocl_ps_sumet = new TH1F("h_hlt_topocl_ps_sumet", "HLT (topocl_PS) Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // HLT SumET [GeV]
  h_hlt_topocl_ps_phi = new TH1F("h_hlt_topocl_ps_phi", "HLT (topocl_PS) MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // HLT phi [GeV]
  wk()->addOutput (h_hlt_topocl_ps_ex);
  wk()->addOutput (h_hlt_topocl_ps_ey);
  wk()->addOutput (h_hlt_topocl_ps_met);
  wk()->addOutput (h_hlt_topocl_ps_sumet);
  wk()->addOutput (h_hlt_topocl_ps_phi);

  // HLT Topocl_puc MET
  h_hlt_topocl_puc_ex = new TH1F("h_hlt_topocl_puc_ex", "HLT (topocl_PUC) Missing E_{x};E_{x} (GeV)", 150, -150,  150); // HLT MEx [GeV]
  h_hlt_topocl_puc_ey = new TH1F("h_hlt_topocl_puc_ey", "HLT (topocl_PUC) Missing E_{y};E_{y} (GeV)", 150, -150,  150); // HLT MEy [GeV]
  h_hlt_topocl_puc_met = new TH1F("h_hlt_topocl_puc_met", "HLT (topocl_PUC) |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // HLT MET [GeV]
  h_hlt_topocl_puc_sumet = new TH1F("h_hlt_topocl_puc_sumet", "HLT (topocl_PUC) Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // HLT SumET [GeV]
  h_hlt_topocl_puc_phi = new TH1F("h_hlt_topocl_puc_phi", "HLT (topocl_PUC) MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // HLT phi [GeV]
  wk()->addOutput (h_hlt_topocl_puc_ex);
  wk()->addOutput (h_hlt_topocl_puc_ey);
  wk()->addOutput (h_hlt_topocl_puc_met);
  wk()->addOutput (h_hlt_topocl_puc_sumet);
  wk()->addOutput (h_hlt_topocl_puc_phi);



  h_jet_selection_pt = new TH1F("h_jet_selection_pt", "Jet Signal p_{T};p_{T} (GeV)", 250, 0, 500); // Jet pt [GeV]
  wk()->addOutput (h_jet_selection_pt);

  h_met_ex = new TH1F("h_met_ex", "Missing E_{x};E_{x} (GeV)", 150, -150,  150); // MEx [GeV]
  h_met_ey = new TH1F("h_met_ey", "Missing E_{y};E_{y} (GeV)", 150, -150,  150); // MEy [GeV]
  h_met = new TH1F("h_met", "|Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // MET [GeV]
  h_sumet = new TH1F("h_sumet", "Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // SumET [GeV]
  h_met_phi = new TH1F("h_met_phi", "MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // MET phi [GeV]
  wk()->addOutput (h_met_ex);
  wk()->addOutput (h_met_ey);
  wk()->addOutput (h_met);
  wk()->addOutput (h_sumet);
  wk()->addOutput (h_met_phi);

  h_emulmet_ex = new TH1F("h_emulmet_ex", "Emulated Missing E_{x};E_{x} (GeV)", 150, -150,  150); // MEx [GeV]
  h_emulmet_ey = new TH1F("h_emulmet_ey", "Emulated Missing E_{y};E_{y} (GeV)", 150, -150,  150); // MEy [GeV]
  h_emulmet = new TH1F("h_emulmet", "Emulated |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // MET [GeV]
  h_emulsumet = new TH1F("h_emulsumet", "Emulated Sum |E_{T}|;SumE_{T} (GeV)", 200, 0, 2000); // SumET [GeV]
  h_emulmet_phi = new TH1F("h_emulmet_phi", "Emulated MET #phi (rad);#phi (rad)", 32, -3.1416, 3.1416); // MET phi [GeV]
  wk()->addOutput (h_emulmet_ex);
  wk()->addOutput (h_emulmet_ey);
  wk()->addOutput (h_emulmet);
  wk()->addOutput (h_emulsumet);
  wk()->addOutput (h_emulmet_phi);


  // Resolution
  /*
  // Define variable bin
  // Increasing bin linearly
  const int nbins = 15;
  double xmin = 0;
  double xmax = 2e3;
  double binwidth = 2*xmax/(nbins*(nbins+1));
  double xbins[nbins+1];
  xbins[0] = xmin;
  for (int i=1;i<=nbins;i++) {
  xbins[i] = xbins[i-1] + binwidth*i;
  }
  */
/*
  // Increasing bin logarithmically
  const int nbins = 15;
  double xmin = 1.;
  double xmax = 2e3;
  double logxmin = log10(xmin);
  double logxmax = log10(xmax);
  double binwidth = (logxmax-logxmin)/nbins;
  double xbins[nbins+1];
  xbins[0] = xmin;
  for (int i=1;i<=nbins;i++) {
    xbins[i] = xmin + pow(10,logxmin+i*binwidth);
  }

  // HLT MEx vs Offline SumET
  h_hlt_ex_offline_sumet = new TH2F("h_hlt_ex_offline_sumet", "HLT (CELL) MET Resolution ;Offline SumE_{T} [GeV];HLT (CELL) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_mht_ex_offline_sumet = new TH2F("h_hlt_mht_ex_offline_sumet", "HLT (MHT) MET Resolution ;Offline SumE_{T} [GeV];HLT (MHT) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_topocl_ex_offline_sumet = new TH2F("h_hlt_topocl_ex_offline_sumet", "HLT (Topocl) MET Resolution ;Offline SumE_{T} [GeV];HLT (Topocl) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_topocl_ps_ex_offline_sumet = new TH2F("h_hlt_topocl_ps_ex_offline_sumet", "HLT (Topocl_ps) MET Resolution ;Offline SumE_{T} [GeV];HLT (Topocl_ps) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_topocl_puc_ex_offline_sumet = new TH2F("h_hlt_topocl_puc_ex_offline_sumet", "HLT (Topocl_puc) MET Resolution ;Offline SumE_{T} [GeV];HLT (Topocl_puc) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  wk()->addOutput (h_hlt_ex_offline_sumet);
  wk()->addOutput (h_hlt_mht_ex_offline_sumet);
  wk()->addOutput (h_hlt_topocl_ex_offline_sumet);
  wk()->addOutput (h_hlt_topocl_ps_ex_offline_sumet);
  wk()->addOutput (h_hlt_topocl_puc_ex_offline_sumet);
  // HLT MEx vs HLT SumET
  h_hlt_ex_hlt_sumet = new TH2F("h_hlt_ex_hlt_sumet", "HLT (CELL) MET Resolution ;HLT (CELL) SumE_{T} [GeV];HLT (CELL) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_mht_ex_hlt_mht_sumet = new TH2F("h_hlt_mht_ex_hlt_mht_sumet", "HLT (MHT) MET Resolution ;HLT (MHT) SumE_{T} [GeV];HLT (MHT) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_topocl_ex_hlt_topocl_sumet = new TH2F("h_hlt_topocl_ex_hlt_topocl_sumet", "HLT (Topocl) MET Resolution ;HLT (Topocl) SumE_{T} [GeV];HLT (Topocl) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_topocl_ps_ex_hlt_topocl_ps_sumet = new TH2F("h_hlt_topocl_ps_ex_hlt_topocl_ps_sumet", "HLT (Topocl_ps) MET Resolution ;HLT (Topocl) SumE_{T} [GeV];HLT (Topocl_ps) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  h_hlt_topocl_puc_ex_hlt_topocl_puc_sumet = new TH2F("h_hlt_topocl_puc_ex_hlt_topocl_puc_sumet", "HLT (Topocl) MET Resolution ;HLT (Topocl_puc) SumE_{T} [GeV];HLT (Topocl_puc) Missing E_{x} [GeV]",nbins,xbins,150,-150,150);
  wk()->addOutput (h_hlt_ex_hlt_sumet);
  wk()->addOutput (h_hlt_mht_ex_hlt_mht_sumet);
  wk()->addOutput (h_hlt_topocl_ex_hlt_topocl_sumet);
  wk()->addOutput (h_hlt_topocl_ps_ex_hlt_topocl_ps_sumet);
  wk()->addOutput (h_hlt_topocl_puc_ex_hlt_topocl_puc_sumet);

  // Linearity
  h_hlt_lin = new TH2F("h_hlt_lin", "HLT (CELL) MET Linearity ;Offline Missing E_{T} [GeV];HLT (CELL) MET / Offline MET", 125, 0, 500, 500, 0, 50);
  h_hlt_mht_lin = new TH2F("h_hlt_mht_lin", "HLT (MHT) MET Linearity ;Offline Missing E_{T} [GeV];HLT (MHT) MET / Offline MET", 125, 0, 500, 500, 0, 50);
  h_hlt_topocl_lin = new TH2F("h_hlt_topocl_lin", "HLT (Topocl) MET Linearity ;Offline Missing E_{T} [GeV];HLT (Topocl) MET / Offline MET", 125, 0, 500, 500, 0, 50);
  h_hlt_topocl_ps_lin = new TH2F("h_hlt_topocl_ps_lin", "HLT (Topocl_ps) MET Linearity ;Offline Missing E_{T} [GeV];HLT (Topocl_ps) MET / Offline MET", 125, 0, 500, 500, 0, 50);
  h_hlt_topocl_puc_lin = new TH2F("h_hlt_topocl_puc_lin", "HLT (Topocl_puc) MET Linearity ;Offline Missing E_{T} [GeV];HLT (Topocl_puc) MET / Offline MET", 125, 0, 500, 500, 0, 50);
  wk()->addOutput (h_hlt_lin);
  wk()->addOutput (h_hlt_mht_lin);
  wk()->addOutput (h_hlt_topocl_lin);
  wk()->addOutput (h_hlt_topocl_ps_lin);
  wk()->addOutput (h_hlt_topocl_puc_lin);

  // HLT Trigger study
  // No selection
  // Calculation of HLT trigger rates
  h_hlt_met_pass_XE50 = new TH1F("h_hlt_met_pass_XE50", "HLT (CELL) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_mht_met_pass_XE50 = new TH1F("h_hlt_mht_met_pass_XE50", "HLT (MHT) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_topocl_met_pass_XE50 = new TH1F("h_hlt_topocl_met_pass_XE50", "HLT (topocl) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_topocl_ps_met_pass_XE50 = new TH1F("h_hlt_topocl_ps_met_pass_XE50", "HLT (topocl_ps) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_topocl_puc_met_pass_XE50 = new TH1F("h_hlt_topocl_puc_met_pass_XE50", "HLT (topocl_puc) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_met_pass_XE60 = new TH1F("h_hlt_met_pass_XE60", "HLT (CELL) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_mht_met_pass_XE60 = new TH1F("h_hlt_mht_met_pass_XE60", "HLT (MHT) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_topocl_met_pass_XE60 = new TH1F("h_hlt_topocl_met_pass_XE60", "HLT (topocl) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_topocl_ps_met_pass_XE60 = new TH1F("h_hlt_topocl_ps_met_pass_XE60", "HLT (topocl_ps) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  h_hlt_topocl_puc_met_pass_XE60 = new TH1F("h_hlt_topocl_puc_met_pass_XE60", "HLT (topocl_puc) Missing E_{x};ME_{T} (GeV)", 250, 0, 500);
  wk()->addOutput (h_hlt_met_pass_XE50);
  wk()->addOutput (h_hlt_mht_met_pass_XE50);
  wk()->addOutput (h_hlt_topocl_met_pass_XE50);
  wk()->addOutput (h_hlt_topocl_ps_met_pass_XE50);
  wk()->addOutput (h_hlt_topocl_puc_met_pass_XE50);
  wk()->addOutput (h_hlt_met_pass_XE60);
  wk()->addOutput (h_hlt_mht_met_pass_XE60);
  wk()->addOutput (h_hlt_topocl_met_pass_XE60);
  wk()->addOutput (h_hlt_topocl_ps_met_pass_XE60);
  wk()->addOutput (h_hlt_topocl_puc_met_pass_XE60);

  // Turn-on Curves
  h_offline_met_pass_hlt_xe70 = new TH1F("h_offline_met_pass_hlt_xe70", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_hlt_xe90 = new TH1F("h_offline_met_pass_hlt_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE50_hlt_xe70 = new TH1F("h_offline_met_pass_l1_XE50_hlt_xe70", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE70_hlt_xe90 = new TH1F("h_offline_met_pass_l1_XE70_hlt_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  wk()->addOutput (h_offline_met_pass_hlt_xe70);
  wk()->addOutput (h_offline_met_pass_hlt_xe90);
  wk()->addOutput (h_offline_met_pass_l1_XE50_hlt_xe70);
  wk()->addOutput (h_offline_met_pass_l1_XE70_hlt_xe90);

  h_offline_met_pass_hlt_mht_xe90 = new TH1F("h_offline_met_pass_hlt_mht_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_hlt_mht_xe110 = new TH1F("h_offline_met_pass_hlt_mht_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE50_hlt_mht_xe90 = new TH1F("h_offline_met_pass_l1_XE50_hlt_mht_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE70_hlt_mht_xe110 = new TH1F("h_offline_met_pass_l1_XE70_hlt_mht_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  wk()->addOutput (h_offline_met_pass_hlt_mht_xe90);
  wk()->addOutput (h_offline_met_pass_hlt_mht_xe110);
  wk()->addOutput (h_offline_met_pass_l1_XE50_hlt_mht_xe90);
  wk()->addOutput (h_offline_met_pass_l1_XE70_hlt_mht_xe110);

  h_offline_met_pass_hlt_topocl_xe90 = new TH1F("h_offline_met_pass_hlt_topocl_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_hlt_topocl_xe110 = new TH1F("h_offline_met_pass_hlt_topocl_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE50_hlt_topocl_xe90 = new TH1F("h_offline_met_pass_l1_XE50_hlt_topocl_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE70_hlt_topocl_xe110 = new TH1F("h_offline_met_pass_l1_XE70_hlt_topocl_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  wk()->addOutput (h_offline_met_pass_hlt_topocl_xe90);
  wk()->addOutput (h_offline_met_pass_hlt_topocl_xe110);
  wk()->addOutput (h_offline_met_pass_l1_XE50_hlt_topocl_xe90);
  wk()->addOutput (h_offline_met_pass_l1_XE70_hlt_topocl_xe110);

  h_offline_met_pass_hlt_topocl_ps_xe90 = new TH1F("h_offline_met_pass_hlt_topocl_ps_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_hlt_topocl_ps_xe110 = new TH1F("h_offline_met_pass_hlt_topocl_ps_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE50_hlt_topocl_ps_xe90 = new TH1F("h_offline_met_pass_l1_XE50_hlt_topocl_ps_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE70_hlt_topocl_ps_xe110 = new TH1F("h_offline_met_pass_l1_XE70_hlt_topocl_ps_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  wk()->addOutput (h_offline_met_pass_hlt_topocl_ps_xe90);
  wk()->addOutput (h_offline_met_pass_hlt_topocl_ps_xe110);
  wk()->addOutput (h_offline_met_pass_l1_XE50_hlt_topocl_ps_xe90);
  wk()->addOutput (h_offline_met_pass_l1_XE70_hlt_topocl_ps_xe110);

  h_offline_met_pass_hlt_topocl_puc_xe90 = new TH1F("h_offline_met_pass_hlt_topocl_puc_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_hlt_topocl_puc_xe110 = new TH1F("h_offline_met_pass_hlt_topocl_puc_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE50_hlt_topocl_puc_xe90 = new TH1F("h_offline_met_pass_l1_XE50_hlt_topocl_puc_xe90", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  h_offline_met_pass_l1_XE70_hlt_topocl_puc_xe110 = new TH1F("h_offline_met_pass_l1_XE70_hlt_topocl_puc_xe110", "Offline |Missing E_{T}|;ME_{T} (GeV)", 250, 0, 500); // Offline MET [GeV]
  wk()->addOutput (h_offline_met_pass_hlt_topocl_puc_xe90);
  wk()->addOutput (h_offline_met_pass_hlt_topocl_puc_xe110);
  wk()->addOutput (h_offline_met_pass_l1_XE50_hlt_topocl_puc_xe90);
  wk()->addOutput (h_offline_met_pass_l1_XE70_hlt_topocl_puc_xe110);

  // Correlation plots
  // L1 vs Offline MET
  h_corr_met_l1_offline = new TH2F("h_corr_met_l1_offline", "L1 vs Offline |Missing E_{T}|;L1 E_{T}^{miss} [GeV];Offline E_{T}^{miss} [GeV]",250,0,500,250,0,500);
  wk()->addOutput (h_corr_met_l1_offline);
  // HLT (CELL) vs Offline MET
  h_corr_met_hlt_offline = new TH2F("h_corr_met_hlt_offline", "HLT (CELL) vs Offline |Missing E_{T}|;HLT (cell) E_{T}^{miss} [GeV];Offline E_{T}^{miss} [GeV]",250,0,500,250,0,500);
  wk()->addOutput (h_corr_met_hlt_offline);
  // HLT (mht) vs Offline MET
  h_corr_met_hlt_mht_offline = new TH2F("h_corr_met_hlt_mht_offline", "HLT (mht) vs Offline |Missing E_{T}|;HLT (mht) E_{T}^{miss} [GeV];Offline E_{T}^{miss} [GeV]",250,0,500,250,0,500);
  wk()->addOutput (h_corr_met_hlt_mht_offline);
  // HLT (topocl) vs Offline MET
  h_corr_met_hlt_topocl_offline = new TH2F("h_corr_met_hlt_topocl_offline", "HLT (topocl) vs Offline |Missing E_{T}|;HLT (topocl) E_{T}^{miss} [GeV];Offline E_{T}^{miss} [GeV]",250,0,500,250,0,500);
  wk()->addOutput (h_corr_met_hlt_topocl_offline);
  // HLT (topocl_ps) vs Offline MET
  h_corr_met_hlt_topocl_ps_offline = new TH2F("h_corr_met_hlt_topocl_ps_offline", "HLT (topocl_ps) vs Offline |Missing E_{T}|;HLT (topocl_ps) E_{T}^{miss} [GeV];Offline E_{T}^{miss} [GeV]",250,0,500,250,0,500);
  wk()->addOutput (h_corr_met_hlt_topocl_ps_offline);
  // HLT (topocl_puc) vs Offline MET
  h_corr_met_hlt_topocl_puc_offline = new TH2F("h_corr_met_hlt_topocl_puc_offline", "HLT (topocl_puc) vs Offline |Missing E_{T}|;HLT (topocl_puc) E_{T}^{miss} [GeV];Offline E_{T}^{miss} [GeV]",250,0,500,250,0,500);
  wk()->addOutput (h_corr_met_hlt_topocl_puc_offline);

*/


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  m_event = wk()->xaodEvent(); // you should have already added this as described before
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  //----------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection in initialise. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // check if the event is data or MC
  // (many tools are applied either to data or MC)
  m_isData = true;
  // check if the event is MC
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    m_isData = false; // can do something with this later
  }

  // count number of events
  m_eventCounter = 0;
  m_numCleanEvents = 0;

  // Enable Cutflow plot
  m_useBitsetCutflow = true;

  // Event Channel
  m_isZvv = true;
  m_isZmumu = false;
  m_isWmunu = false;
  m_isZee = false;
  m_isWenu = false;

  // Enable Overlap Removal tool
  m_doORtool = false;
  m_doORmanual = true;

  // Cut values
  m_muonPtCut = 7.; /// GeV
  m_muonEtaCut = 2.5;
  m_elecPtCut = 7.; /// GeV
  m_elecEtaCut = 2.47;
  m_photPtCut = 20.; /// GeV
  m_photEtaCut = 2.47;
  m_jetPtCut = 20.; /// GeV
  m_jetEtaCut = 4.5;
  m_diJet1PtCut = 80.; /// GeV
  m_diJet2PtCut = 50.; /// GeV
  m_diJetEtaCut = 4.4;
  m_CJVptCut = 25.; ///GeV
  m_metCut = 200.; ///GeV
  m_mjjCut = 200.; ///GeV
  m_LeadLepPtCut = 80.; ///GeV
  m_SubLeadLepPtCut = 7.; ///GeV
  m_ORJETdeltaR = 0.2;
  m_isoMuonPtMin = 10.; ///GeV
  m_isoMuonPtMax = 500.; ///GeV
  m_recoSF = true;
  m_isoSF = true;
  m_ttvaSF = true;

  // GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  std::vector<std::string> vecStringGRL;
  // GRL xml file should be put in MetTriggerPackage/share directory
  vecStringGRL.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/MetTriggerPackage/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml"));
  EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
  EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
  EL_RETURN_CHECK("initialize()",m_grl->initialize());

  // Initialize and configure trigger tools
  m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
  EL_RETURN_CHECK( "initialize", m_trigConfigTool->initialize() );
  ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
  m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ) ); // connect the TrigDecisionTool to the ConfigTool
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision" ) );
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->initialize() );

  // initialize the muon calibration and smearing tool
  m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" );
  //m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::DEBUG );
  m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::INFO );
  EL_RETURN_CHECK("initialize()",m_muonCalibrationAndSmearingTool->initialize());

  // initialize the electron and photon calibration and smearing tool
  m_egammaCalibrationAndSmearingTool = new CP::EgammaCalibrationAndSmearingTool( "EgammaCorrectionTool" );
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "ESModel", "es2015PRE" ));  // see below for options
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "decorrelationModel", "FULL_ETACORRELATED_v1" ));  // see below for options
  //EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->setProperty( "decorrelationModel", "FULL_v1" ));  // see below for options
  EL_RETURN_CHECK("initialize()",m_egammaCalibrationAndSmearingTool->initialize());

  // Initialize the MC fudge tool
  m_electronPhotonShowerShapeFudgeTool = new ElectronPhotonShowerShapeFudgeTool( "ElectronPhotonShowerShapeFudgeTool" );
  int FFset = 16; // for MC15 samples, which are based on a geometry derived from GEO-21
  EL_RETURN_CHECK("initialize()",m_electronPhotonShowerShapeFudgeTool->setProperty("Preselection", FFset));
  EL_RETURN_CHECK("initialize()",m_electronPhotonShowerShapeFudgeTool->initialize() );

  // Muon identification (Medium)
  // initialize the muon selection tool
  m_muonSelection = new CP::MuonSelectionTool( "MuonSelection" );
  //EL_RETURN_CHECK("initialize()",m_muonSelection->setProperty( "MaxEta", 2.5 ));
  EL_RETURN_CHECK("initialize()",m_muonSelection->setProperty( "MuQuality", 1)); // 0 tight, 1 medium, 2 loose, 3 very loose
  //m_muonSelection->msg().setLevel( MSG::VERBOSE );
  m_muonSelection->msg().setLevel( MSG::INFO );
  //m_muonSelection->msg().setLevel( MSG::ERROR );
  EL_RETURN_CHECK("initialize()",m_muonSelection->initialize());
  // Muon identification (Loose)
  m_loosemuonSelection = new CP::MuonSelectionTool( "MuonLooseSelection" );
  //m_loosemuonSelection->msg().setLevel( MSG::VERBOSE );
  m_loosemuonSelection->msg().setLevel( MSG::INFO );
  //m_loosemuonSelection->msg().setLevel( MSG::ERROR );
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->setProperty( "MaxEta", 2.5 ));
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->setProperty( "MuQuality", 2));
  EL_RETURN_CHECK("initialize()",m_loosemuonSelection->initialize());

  // LH Electron identification
  // initialize the electron selection tool
  m_LHToolTight2015    = new AsgElectronLikelihoodTool ("m_LHToolTight2015");
  m_LHToolMedium2015   = new AsgElectronLikelihoodTool ("m_LHToolMedium2015"); 
  m_LHToolLoose2015    = new AsgElectronLikelihoodTool ("m_LHToolLoose2015");
  // initialize the primary vertex container for the tool to have access to the number of vertices used to adapt cuts based on the pileup
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  // define the config files
  std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150712/";
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf"));
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
  //EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2015.conf"));
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2015_CutBL.conf"));
  // initialize
  EL_RETURN_CHECK("initialize()",m_LHToolTight2015->initialize());
  EL_RETURN_CHECK("initialize()",m_LHToolMedium2015->initialize());
  EL_RETURN_CHECK("initialize()",m_LHToolLoose2015->initialize());

  // Photon identification (Tight)
  // create the selector
  m_photonTightIsEMSelector = new AsgPhotonIsEMSelector ( "PhotonTightIsEMSelector" );
  // decide which kind of selection (Loose/Medium/Tight) you want to use
  EL_RETURN_CHECK("initialize()",m_photonTightIsEMSelector->setProperty("isEMMask",egammaPID::PhotonTight));
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  EL_RETURN_CHECK("initialize()",m_photonTightIsEMSelector->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf"));
  // initialise the tool
  EL_RETURN_CHECK("initialize()",m_photonTightIsEMSelector->initialize());
  // Photon identification (Medium)
  // create the selector
  m_photonMediumIsEMSelector = new AsgPhotonIsEMSelector ( "PhotonMediumIsEMSelector" );
  // decide which kind of selection (Loose/Medium/Tight) you want to use
  EL_RETURN_CHECK("initialize()",m_photonMediumIsEMSelector->setProperty("isEMMask",egammaPID::PhotonMedium));
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  EL_RETURN_CHECK("initialize()",m_photonMediumIsEMSelector->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMMediumSelectorCutDefs.conf"));
  // initialise the tool
  EL_RETURN_CHECK("initialize()",m_photonMediumIsEMSelector->initialize());
  // Photon identification (Loose)
  // create the selector
  m_photonLooseIsEMSelector = new AsgPhotonIsEMSelector ( "PhotonLooseIsEMSelector" );
  // decide which kind of selection (Loose/Medium/Tight) you want to use
  EL_RETURN_CHECK("initialize()",m_photonLooseIsEMSelector->setProperty("isEMMask",egammaPID::PhotonLoose));
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  EL_RETURN_CHECK("initialize()",m_photonLooseIsEMSelector->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMLooseSelectorCutDefs.conf"));
  // initialise the tool
  EL_RETURN_CHECK("initialize()",m_photonLooseIsEMSelector->initialize());



  // IsolationSelectionTool
  m_IsolationSelectionTool = new CP::IsolationSelectionTool("IsolationSelectionTool");
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->setProperty("MuonWP","Gradient"));
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->setProperty("ElectronWP","Gradient"));
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->setProperty("PhotonWP","Cone40"));
  EL_RETURN_CHECK("initialize()",m_IsolationSelectionTool->initialize());
  // IsolationSelectionTool for VBF signal
  m_IsoToolVBF = new CP::IsolationSelectionTool("IsoToolVBF");
  //EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("MuonWP","FixedCutLoose"));
  //EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("ElectronWP","FixedCutLoose"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("MuonWP","LooseTrackOnly"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("ElectronWP","LooseTrackOnly"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->setProperty("PhotonWP","Cone40"));
  EL_RETURN_CHECK("initialize()",m_IsoToolVBF->initialize());


  // Tau identification
  // initialize the tau selection tool
  m_tauSelTool = new TauAnalysisTools::TauSelectionTool( "TauSelectionTool" );
  // define the config files
  std::string confPath = "$ROOTCOREBIN/data/MetTriggerPackage/";
  EL_RETURN_CHECK("initialize()",m_tauSelTool->setProperty( "ConfigPath", confPath+"recommended_selection_mc15.conf"));
  m_tauSelTool->msg().setLevel( MSG::INFO );
  //m_tauSelTool->msg().setLevel( MSG::DEBUG );
  // initialize
  EL_RETURN_CHECK("initialize()",m_tauSelTool->initialize());

  // initialize the tau selection tool for VBF analysis
  m_tauSelToolVBF = new TauAnalysisTools::TauSelectionTool( "TauSelectionToolVBF" );
  // define the config files
  EL_RETURN_CHECK("initialize()",m_tauSelToolVBF->setProperty( "ConfigPath", confPath+"recommended_selection_mc15_VBF.conf"));
  m_tauSelToolVBF->msg().setLevel( MSG::INFO );
  // initialize
  EL_RETURN_CHECK("initialize()",m_tauSelToolVBF->initialize());

  // Initialise tau smearing tool
  m_tauSmearingTool = new TauAnalysisTools::TauSmearingTool( "TauSmaringTool" );
  m_tauSmearingTool->msg().setLevel( MSG::INFO );
  // initialize
  EL_RETURN_CHECK("initialize()",m_tauSmearingTool->initialize());

  // Initialise TauOverlappingElectronLLHDecorator
  m_tauOverlappingElectronLLHDecorator = new TauAnalysisTools::TauOverlappingElectronLLHDecorator("TauOverlappingElectronLLHDecorator"); 
  EL_RETURN_CHECK("initialize()",m_tauOverlappingElectronLLHDecorator->initialize());


  // Jet
  // JES Calibration (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JES_calibration_AN1)
  const std::string name = "MetTrigxAODAnalysis"; //string describing the current thread, for logging
  TString jetAlgo = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo
  //TString config = "JES_MC15Prerecommendation_April2015.config"; //Path to global config used to initialize the tool
  TString config = "JES_2015dataset_recommendation_Feb2016.config"; //Path to global config used to initialize the tool
  TString calibSeq = "JetArea_Residual_Origin_EtaJES_GSC"; //String describing the calibration sequence to apply
  if (m_isData) calibSeq += "_Insitu";
  //Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
  m_jetCalibration = new JetCalibrationTool(name, jetAlgo, config, calibSeq, m_isData);
  //Initialize the tool
  EL_RETURN_CHECK("initialize()",m_jetCalibration->initializeTool(name));

  // JES uncertainty (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JES_uncertainty)
  m_jetUncertaintiesTool = new JetUncertaintiesTool("JetUncertaintiesTool");
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("JetDefinition", "AntiKt4EMTopo"));
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("MCType", "MC15"));
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->setProperty("ConfigFile", "JES_2015/Prerec/PrerecJES2015_AllNuisanceParameters_25ns.config"));
  // Initialise jet uncertainty tool
  EL_RETURN_CHECK("initialize()",m_jetUncertaintiesTool->initialize());

  // JER uncertainty  (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15#JER_uncertainty)
  // Jet Resolution (https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetResolution2015Prerecom)
  // Configure the JERTool
  m_jerTool = new JERTool("JERTool");
  EL_RETURN_CHECK("initialize()",m_jerTool->setProperty("PlotFileName", "JetResolution/Prerec2015_xCalib_2012JER_ReducedTo9NP_Plots_v2.root"));
  EL_RETURN_CHECK("initialize()",m_jerTool->setProperty("CollectionName", "AntiKt4EMTopoJets"));
  EL_RETURN_CHECK("initialize()",m_jerTool->initialize());
  // Configure the JERSmearingTool
  m_jerHandle = ToolHandle<IJERTool>(m_jerTool->name());
  m_jerSmearingTool = new JERSmearingTool("JERSmearingTool");
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("ApplyNominalSmearing", false));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("JERTool", m_jerHandle));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("isMC", !m_isData));
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->setProperty("SystematicMode", "Simple")); //"Simple" provides one NP (smearing only in MC), "Full" provides 10NPs (smearing both on data and MC)
  EL_RETURN_CHECK("initialize()",m_jerSmearingTool->initialize());

  // Configure the JVT tool.
  m_jvtag = 0;
  m_jvtag = new JetVertexTaggerTool("jvtag");
  m_jvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
  EL_RETURN_CHECK("initialize()",m_jvtag->setProperty("JVTFileName","JetMomentTools/JVTlikelihood_20140805.root"));
  EL_RETURN_CHECK("initialize()",m_jvtag->initialize());

  // Initialize and configure the jet cleaning tool
  m_jetCleaningTight = new JetCleaningTool("JetCleaningTight");
  m_jetCleaningLoose = new JetCleaningTool("JetCleaningLoose");
  m_jetCleaningTight->msg().setLevel( MSG::DEBUG ); 
  m_jetCleaningLoose->msg().setLevel( MSG::DEBUG ); 
  EL_RETURN_CHECK("initialize()",m_jetCleaningTight->setProperty( "CutLevel", "TightBad"));
  EL_RETURN_CHECK("initialize()",m_jetCleaningLoose->setProperty( "CutLevel", "LooseBad"));
  //EL_RETURN_CHECK("initialize()",m_jetCleaningTight->setProperty("DoUgly", false));
  //EL_RETURN_CHECK("initialize()",m_jetCleaningLoose->setProperty("DoUgly", false));
  EL_RETURN_CHECK("initialize()",m_jetCleaningTight->initialize());
  EL_RETURN_CHECK("initialize()",m_jetCleaningLoose->initialize());

  // Initialise MET tools
  m_metMaker = new met::METMaker("METMakerTool");
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoRemoveMuonJets", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoSetMuonJetEMScale", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoMuonEloss", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty( "DoIsolMuonEloss", true));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty("JetMinWeightedPt", 20000.));
  //EL_RETURN_CHECK("initialize()",m_metMaker->setProperty("JetMinEFrac", 0.0));
  //m_metMaker->msg().setLevel( MSG::VERBOSE ); // or DEBUG or VERBOSE
  EL_RETURN_CHECK("initialize()",m_metMaker->initialize());

  // Initialize the harmonization reccommendation tools
  const bool doTaus = true, doPhotons = true;
  const bool boostedLeptons = false;
  EL_RETURN_CHECK("initialize()",ORUtils::recommendedTools(m_toolBox, "OverlapRemovalTool", 
                                                          inputLabel, outputLabel, bJetLabel, 
                                                          boostedLeptons, outputPassValue, 
                                                          doTaus, doPhotons));
  //EL_RETURN_CHECK("initialize()",ORUtils::harmonizedTools(m_toolBox, "OverlapRemovalTool", 
  //      inputLabel, outputLabel,
  //      outputPassValue, doTaus, doPhotons));
  // Set message level for all tools
  //m_toolBox->setMsgLevel(MSG::DEBUG);
  // Initialize all tools
  m_orTool = static_cast<ORUtils::OverlapRemovalTool*>(m_toolBox.getMasterTool());
  m_orTool->setName("ORTool");
  EL_RETURN_CHECK("initialize()",m_toolBox.initialize());


  // Initialise Muon Efficiency Tool
  m_muonEfficiencySFTool = new CP::MuonEfficiencyScaleFactors( "MuonEfficiencySFTool" );
  EL_RETURN_CHECK("initialize()",m_muonEfficiencySFTool->setProperty("WorkingPoint", "Loose") );
  EL_RETURN_CHECK("initialize()",m_muonEfficiencySFTool->setProperty("CalibrationRelease", "Data15_allPeriods_260116"));
  EL_RETURN_CHECK("initialize()",m_muonEfficiencySFTool->initialize() );
  // Initialise Muon Isolation Tool
  m_muonIsolationSFTool = new CP::MuonEfficiencyScaleFactors( "MuonIsolationSFTool" );
  EL_RETURN_CHECK("initialize()",m_muonIsolationSFTool->setProperty("WorkingPoint", "LooseTrackOnlyIso") );
  EL_RETURN_CHECK("initialize()",m_muonIsolationSFTool->setProperty("CalibrationRelease", "Data15_allPeriods_260116"));
  EL_RETURN_CHECK("initialize()",m_muonIsolationSFTool->initialize() );
  // Initialise Muon TTVA Efficiency Tool
  m_muonTTVAEfficiencySFTool = new CP::MuonEfficiencyScaleFactors( "MuonTTVAEfficiencySFTool" );
  EL_RETURN_CHECK("initialize()",m_muonTTVAEfficiencySFTool->setProperty("WorkingPoint", "TTVA") );
  EL_RETURN_CHECK("initialize()",m_muonTTVAEfficiencySFTool->setProperty("CalibrationRelease", "Data15_allPeriods_260116"));
  EL_RETURN_CHECK("initialize()",m_muonTTVAEfficiencySFTool->initialize() );

  // Initialise Muon Trigger Scale Factor Tool
  m_muonTriggerSFTool = new CP::MuonTriggerScaleFactors( "MuonTriggerSFTool" );
  EL_RETURN_CHECK("initialize()",m_muonTriggerSFTool->setProperty("MuonQuality", "Loose"));
  EL_RETURN_CHECK("initialize()",m_muonTriggerSFTool->setProperty("Isolation", "LooseTrackOnly"));
  EL_RETURN_CHECK("initialize()",m_muonTriggerSFTool->initialize() );

  // Initialise Electron Efficiency Tool
  m_elecEfficiencySFTool_reco = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_reco");
  std::vector< std::string > corrFileNameList_reco;
  corrFileNameList_reco.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.RecoTrk.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("CorrectionFileNameList", corrFileNameList_reco) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_reco->initialize() );

  m_elecEfficiencySFTool_id_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Loose");
  std::vector< std::string > corrFileNameList_id_Loose;
  corrFileNameList_id_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.LooseAndBLayerLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("CorrectionFileNameList", corrFileNameList_id_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Loose->initialize() );

  m_elecEfficiencySFTool_id_Medium = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Medium");
  std::vector< std::string > corrFileNameList_id_Medium;
  corrFileNameList_id_Medium.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.MediumLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("CorrectionFileNameList", corrFileNameList_id_Medium) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Medium->initialize() );

  m_elecEfficiencySFTool_id_Tight = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_id_Tight");
  std::vector< std::string > corrFileNameList_id_Tight;
  corrFileNameList_id_Tight.push_back("ElectronEfficiencyCorrection/efficiencySF.offline.TightLLH_d0z0.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->setProperty("CorrectionFileNameList", corrFileNameList_id_Tight) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_id_Tight->initialize() );

  m_elecEfficiencySFTool_iso_Loose = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Loose");
  std::vector< std::string > corrFileNameList_iso_Loose;
  corrFileNameList_iso_Loose.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("CorrectionFileNameList", corrFileNameList_iso_Loose) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Loose->initialize() );

  m_elecEfficiencySFTool_iso_Medium = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Medium");
  std::vector< std::string > corrFileNameList_iso_Medium;
  corrFileNameList_iso_Medium.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.MediumLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("CorrectionFileNameList", corrFileNameList_iso_Medium) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->setProperty("ForceDataType", 1) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Medium->initialize() );

  m_elecEfficiencySFTool_iso_Tight = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_iso_Tight");
  std::vector< std::string > corrFileNameList_iso_Tight;
  corrFileNameList_iso_Tight.push_back("ElectronEfficiencyCorrection/efficiencySF.Isolation.TightLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->setProperty("CorrectionFileNameList", corrFileNameList_iso_Tight) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_iso_Tight->initialize() );

  m_elecEfficiencySFTool_trigEff = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_trigEff");
  std::vector< std::string > corrFileNameList_trigEff;
  corrFileNameList_trigEff.push_back("ElectronEfficiencyCorrection/efficiency.e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->setProperty("CorrectionFileNameList", corrFileNameList_trigEff) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigEff->initialize() );

  m_elecEfficiencySFTool_trigSF = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyCorrectionTool_trigSF");
  std::vector< std::string > corrFileNameList_trigSF;
  corrFileNameList_trigSF.push_back("ElectronEfficiencyCorrection/efficiencySF.e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose.LooseAndBLayerLLH_d0z0_v8_isolLooseTrackOnly.2015.13TeV.rel20p0.25ns.v04.root");
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF->setProperty("CorrectionFileNameList", corrFileNameList_trigSF) );
  EL_RETURN_CHECK("initialize()",m_elecEfficiencySFTool_trigSF->initialize() );

  // Initialise Jet JVT Efficiency Tool
  m_jvtefficiencyTool = new CP::JetJvtEfficiency("JvtEfficiencyTool");
  //EL_RETURN_CHECK("initialize()",m_jvtefficiencyTool->setProperty("WorkingPoint",) );
  EL_RETURN_CHECK("initialize()",m_jvtefficiencyTool->initialize() );

  // Initialise Tau Efficiency Tool
  m_tauEffTool = new TauAnalysisTools::TauEfficiencyCorrectionsTool("TauEffTool");
  EL_RETURN_CHECK("initialize()",m_tauEffTool->initialize() );

  // Initialise MET Tools
  m_metSystTool = new met::METSystematicsTool("METSystTool");
  EL_RETURN_CHECK("initialize()",m_metSystTool->setProperty("JetColl", "AntiKt4EMTopoJets") );
  EL_RETURN_CHECK("initialize()",m_metSystTool->setProperty("ConfigSoftTrkFile", "TrackSoftTerms.config") );
  EL_RETURN_CHECK("initialize()",m_metSystTool->initialize() );

  // Initialise Isolation Correction Tool
  m_isoCorrTool = new CP::IsolationCorrectionTool( "IsoCorrTool" );
  EL_RETURN_CHECK("initialize()",m_isoCorrTool->setProperty( "IsMC", !m_isData) );
  //EL_RETURN_CHECK("initialize()",m_isoCorrTool->setProperty( "AFII_corr", m_isAtlfast) );
  EL_RETURN_CHECK("initialize()",m_isoCorrTool->initialize() );

  // Initialise PileupReweighting Tool
  m_prwTool = new CP::PileupReweightingTool("PrwTool");
  std::vector<std::string> file_conf;
  // xml file should be put in MetTriggerPackage/share directory
  file_conf.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/MetTriggerPackage/PRW.root"));
  std::vector<std::string> file_ilumi;
  file_ilumi.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/MetTriggerPackage/ilumicalc_histograms_None_276262-284484.root"));
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("ConfigFiles", file_conf) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("LumiCalcFiles", file_ilumi) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactor",     1. / 1.16) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactorUP",   1.) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("DataScaleFactorDOWN", 1. / 1.23) );
  EL_RETURN_CHECK("initialize()",m_prwTool->setProperty("UnrepresentedDataAction", 2));
  EL_RETURN_CHECK("initialize()",m_prwTool->initialize() );


  
  // Initialize Cutflow
  if (m_useBitsetCutflow)
    m_BitsetCutflow = new BitsetCutflow(wk());


  // Initialize Cutflow count array
  for (int i=0; i<20; i++) {
    m_eventCutflow[i]=0;
  }  

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // push cutflow bitset to cutflow hist
  if (m_useBitsetCutflow)
    m_BitsetCutflow->PushBitSet();

  // print every 100 events, so we know where we are:
  if( (m_eventCounter % 100) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;
  m_eventCutflow[0]+=1;

  //----------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection in execute. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Calculate EventWeight
  if (m_isData) // it's data!
    mcEventWeight = 1.0;
  else {
    float pu_weight = m_prwTool->getCombinedWeight(*eventInfo); // Get Pile-up weight
    mcEventWeight = eventInfo->mcEventWeight() * pu_weight;
  }

  // MC Channel Number
  if (m_isData)
    mcChannelNumber = 1;
  else mcChannelNumber = eventInfo->mcChannelNumber();

  // BCID Information
  int m_Bcid = 0;
  bool m_passBcid = true;
  if (m_isData){
    m_Bcid = eventInfo->bcid();
    h_bcid->Fill(m_Bcid);
    if ( (m_Bcid == 1) ||
        (m_Bcid >= 170 && m_Bcid <= 181) || (m_Bcid >= 219 && m_Bcid <= 230) ||
        (m_Bcid >= 328 && m_Bcid <= 339) || (m_Bcid >= 437 && m_Bcid <= 448) ||
        (m_Bcid >= 1113 && m_Bcid <= 1124) || (m_Bcid >= 1222 && m_Bcid <= 1233) ||
        (m_Bcid >= 1331 && m_Bcid <= 1342) || (m_Bcid >= 2007 && m_Bcid <= 2018) ||
        (m_Bcid >= 2116 && m_Bcid <= 2127) || (m_Bcid >= 2225 && m_Bcid <= 2236) ||
        (m_Bcid >= 2901 && m_Bcid <= 2912) || (m_Bcid >= 3010 && m_Bcid <= 3021) ||
        (m_Bcid >= 3119 && m_Bcid <= 3130) ){
      m_passBcid = false;
    }
  }


  /*
  // if data check if event passes GRL
  if(m_isData){ // it's data!
    if(!m_grl->passRunLB(*eventInfo)){
      return EL::StatusCode::SUCCESS; // go to next event
    }
  } // end if not MC
  */
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("GRL");
  m_eventCutflow[1]+=1;


  //------------------------------------------------------------
  // Apply event cleaning to remove events due to 
  // problematic regions of the detector, and incomplete events.
  // Apply to data.
  //------------------------------------------------------------
  // reject event if:
  if(m_isData){
    if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) 
        || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) 
        || (eventInfo->errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error) 
        || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) )  )
    {
      return EL::StatusCode::SUCCESS; // go to the next event
    } // end if event flags check
  } // end if the event is data
  m_numCleanEvents++;
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("LAr_Tile_Core");
  m_eventCutflow[2]+=1;


  //---------------------
  // EVENT SELECTION
  //---------------------

  //---------------------
  // Retrive vertex object and select events with at least one good primary vertex with at least 2 tracks
  //---------------------
  const xAOD::VertexContainer* vertices(0);
  /// retrieve arguments: container type, container key
  if ( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){ 
    Error("execute()","Failed to retrieve PrimaryVertices container. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  xAOD::Vertex* primVertex = 0;
  int nGoodVtx = 0;
  for( ; vtx_itr != vtx_end; ++vtx_itr ) {
    if ((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx){
      primVertex = (*vtx_itr);
      nGoodVtx++;
    }
  }


  //--------------
  // Preseletion
  //--------------
  //if (nGoodVtx==1)
  //   Info("execute()", "  Found one prim.vertex: nGoodVtx = %d", nGoodVtx); // just to print out something
  if (nGoodVtx>1)
    Info("execute()", "  WARNING!!!! Found more than one prim.vertex: nGoodVtx = %d", nGoodVtx); // just to print out something
  if (nGoodVtx==0)
    //   Info("execute()", "  %s", "No one prim.vertex found"); // just to print out something
    return EL::StatusCode::SUCCESS;
  if (primVertex->nTrackParticles() < 2) return EL::StatusCode::SUCCESS;

  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Primary vertex");
  m_eventCutflow[3]+=1;


  // Number of Interactions
  float m_AverageInteractionsPerCrossing = eventInfo->averageInteractionsPerCrossing();
  float m_ActualInteractionsPerCrossing = eventInfo->actualInteractionsPerCrossing();
  h_avg_interaction->Fill(m_AverageInteractionsPerCrossing);
  h_act_interaction->Fill(m_ActualInteractionsPerCrossing);


  //-------------------------------- 
  // Retrieve MET Trigger Containers
  //-------------------------------- 

  // get xAOD container of interest

  // retrieve xAOD L1 ROI
  const xAOD::EnergySumRoI *m_l1_roi_cont = nullptr;
  if( ! m_event->retrieve( m_l1_roi_cont, "LVL1EnergySumRoI" ).isSuccess() ){
    Error("execute()", "Failed to retrieve LVL1EnergySumRoI container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // retrieve HLT containers
  // Get HLT (CELL) container
  const xAOD::TrigMissingETContainer *m_hlt_met_cont = nullptr;
  if( ! m_event->retrieve( m_hlt_met_cont, "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET" ).isSuccess() ){
    Error("execute()", "Failed to retrieve TrigEFMissingET container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Get HLT mht container
  const xAOD::TrigMissingETContainer *m_hlt_met_cont_mht = nullptr;
  if( ! m_event->retrieve( m_hlt_met_cont_mht, "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET_mht" ).isSuccess() ){
    Error("execute()", "Failed to retrieve TrigEFMissingET_mht container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Get HLT topocl container
  const xAOD::TrigMissingETContainer *m_hlt_met_cont_topocl = nullptr;
  if( ! m_event->retrieve( m_hlt_met_cont_topocl, "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET_topocl" ).isSuccess() ){
    Error("execute()", "Failed to retrieve TrigEFMissingET_topocl container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Get HLT topocl_PS container
  const xAOD::TrigMissingETContainer *m_hlt_met_cont_topocl_ps = nullptr;
  if( ! m_event->retrieve( m_hlt_met_cont_topocl_ps, "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET_topocl_PS" ).isSuccess() ){
    Error("execute()", "Failed to retrieve TrigEFMissingET_topocl_PS container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // Get HLT topocl_PUC container
  const xAOD::TrigMissingETContainer *m_hlt_met_cont_topocl_puc = nullptr;
  if( ! m_event->retrieve( m_hlt_met_cont_topocl_puc, "HLT_xAOD__TrigMissingETContainer_TrigEFMissingET_topocl_PUC" ).isSuccess() ){
    Error("execute()", "Failed to retrieve TrigEFMissingET_topocl_PUC container. Exiting." );
    return EL::StatusCode::FAILURE;
  }


  // L1
  float l1_mex = -9e9;
  float l1_mey = -9e9;
  float l1_met = -9e9;
  float l1_sumet = -9e9;
  float l1_phi = -9e9;

  l1_mex = - (m_l1_roi_cont->energyX()) * 0.001;
  l1_mey = - (m_l1_roi_cont->energyY()) * 0.001;
  l1_met = sqrt(l1_mex*l1_mex + l1_mey*l1_mey);
  l1_sumet = (m_l1_roi_cont->energyT()) * 0.001;
  l1_phi = atan2f(l1_mey,l1_mex);

  h_l1_mex->Fill( l1_mex ); // GeV
  h_l1_mey->Fill( l1_mey ); // GeV
  h_l1_met->Fill( l1_met ); // GeV
  h_l1_sumet->Fill( l1_sumet ); // GeV
  h_l1_phi->Fill( l1_phi ); // GeV


  // HLT (CELL)
  float hlt_ex = -9e9;
  float hlt_ey = -9e9;
  float hlt_met = -9e9;
  float hlt_sumet = -9e9;
  float hlt_phi = -9e9;

  // Get MET values
  const xAOD::TrigMissingET *m_hlt_met = nullptr;

  m_hlt_met = m_hlt_met_cont->at(0);

  hlt_ex = m_hlt_met->ex() * 0.001; 
  hlt_ey = m_hlt_met->ey() * 0.001; 
  hlt_met = sqrt(hlt_ex*hlt_ex+hlt_ey*hlt_ey);
  hlt_sumet = m_hlt_met->sumEt() * 0.001;
  hlt_phi = atan2(hlt_ey, hlt_ex);
  //   Info("execute()", "  HLT_MET ex (test) = %.2f GeV", hlt_ex); // just to print out something

  h_hlt_ex->Fill( hlt_ex ); // GeV
  h_hlt_ey->Fill( hlt_ey ); // GeV
  h_hlt_met->Fill( hlt_met ); // GeV
  h_hlt_sumet->Fill( hlt_sumet ); // GeV
  h_hlt_phi->Fill( hlt_phi ); // GeV


  // HLT (MHT)
  float hlt_mht_ex = -9e9;
  float hlt_mht_ey = -9e9;
  float hlt_mht_met = -9e9;
  float hlt_mht_sumet = -9e9;
  float hlt_mht_phi = -9e9;

  // Get MET values
  const xAOD::TrigMissingET *m_hlt_met_mht = nullptr;

  m_hlt_met_mht = m_hlt_met_cont_mht->at(0);

  hlt_mht_ex = m_hlt_met_mht->ex() * 0.001; 
  hlt_mht_ey = m_hlt_met_mht->ey() * 0.001; 
  hlt_mht_met = sqrt(hlt_mht_ex*hlt_mht_ex+hlt_mht_ey*hlt_mht_ey);
  hlt_mht_sumet = m_hlt_met_mht->sumEt() * 0.001;
  hlt_mht_phi = atan2(hlt_mht_ey, hlt_mht_ex);

  h_hlt_mht_ex->Fill( hlt_mht_ex ); // GeV
  h_hlt_mht_ey->Fill( hlt_mht_ey ); // GeV
  h_hlt_mht_met->Fill( hlt_mht_met ); // GeV
  h_hlt_mht_sumet->Fill( hlt_mht_sumet ); // GeV
  h_hlt_mht_phi->Fill( hlt_mht_phi ); // GeV


  // HLT (topocl)
  float hlt_topocl_ex = -9e9;
  float hlt_topocl_ey = -9e9;
  float hlt_topocl_met = -9e9;
  float hlt_topocl_sumet = -9e9;
  float hlt_topocl_phi = -9e9;

  // Get MET values
  const xAOD::TrigMissingET *m_hlt_met_topocl = nullptr;

  m_hlt_met_topocl = m_hlt_met_cont_topocl->at(0);

  hlt_topocl_ex = m_hlt_met_topocl->ex() * 0.001; 
  hlt_topocl_ey = m_hlt_met_topocl->ey() * 0.001; 
  hlt_topocl_met = sqrt(hlt_topocl_ex*hlt_topocl_ex+hlt_topocl_ey*hlt_topocl_ey);
  hlt_topocl_sumet = m_hlt_met_topocl->sumEt() * 0.001;
  hlt_topocl_phi = atan2(hlt_topocl_ey, hlt_topocl_ex);

  h_hlt_topocl_ex->Fill( hlt_topocl_ex ); // GeV
  h_hlt_topocl_ey->Fill( hlt_topocl_ey ); // GeV
  h_hlt_topocl_met->Fill( hlt_topocl_met ); // GeV
  h_hlt_topocl_sumet->Fill( hlt_topocl_sumet ); // GeV
  h_hlt_topocl_phi->Fill( hlt_topocl_phi ); // GeV


  // HLT (topocl_PS)
  float hlt_topocl_ps_ex = -9e9;
  float hlt_topocl_ps_ey = -9e9;
  float hlt_topocl_ps_met = -9e9;
  float hlt_topocl_ps_sumet = -9e9;
  float hlt_topocl_ps_phi = -9e9;

  // Get MET values
  const xAOD::TrigMissingET *m_hlt_met_topocl_ps = nullptr;

  m_hlt_met_topocl_ps = m_hlt_met_cont_topocl_ps->at(0);

  hlt_topocl_ps_ex = m_hlt_met_topocl_ps->ex() * 0.001; 
  hlt_topocl_ps_ey = m_hlt_met_topocl_ps->ey() * 0.001; 
  hlt_topocl_ps_met = sqrt(hlt_topocl_ps_ex*hlt_topocl_ps_ex+hlt_topocl_ps_ey*hlt_topocl_ps_ey);
  hlt_topocl_ps_sumet = m_hlt_met_topocl_ps->sumEt() * 0.001;
  hlt_topocl_ps_phi = atan2(hlt_topocl_ps_ey, hlt_topocl_ps_ex);

  h_hlt_topocl_ps_ex->Fill( hlt_topocl_ps_ex ); // GeV
  h_hlt_topocl_ps_ey->Fill( hlt_topocl_ps_ey ); // GeV
  h_hlt_topocl_ps_met->Fill( hlt_topocl_ps_met ); // GeV
  h_hlt_topocl_ps_sumet->Fill( hlt_topocl_ps_sumet ); // GeV
  h_hlt_topocl_ps_phi->Fill( hlt_topocl_ps_phi ); // GeV


  // HLT (topocl_PUC)
  float hlt_topocl_puc_ex = -9e9;
  float hlt_topocl_puc_ey = -9e9;
  float hlt_topocl_puc_met = -9e9;
  float hlt_topocl_puc_sumet = -9e9;
  float hlt_topocl_puc_phi = -9e9;

  // Get MET values
  const xAOD::TrigMissingET *m_hlt_met_topocl_puc = nullptr;

  m_hlt_met_topocl_puc = m_hlt_met_cont_topocl_puc->at(0);

  hlt_topocl_puc_ex = m_hlt_met_topocl_puc->ex() * 0.001; 
  hlt_topocl_puc_ey = m_hlt_met_topocl_puc->ey() * 0.001; 
  hlt_topocl_puc_met = sqrt(hlt_topocl_puc_ex*hlt_topocl_puc_ex+hlt_topocl_puc_ey*hlt_topocl_puc_ey);
  hlt_topocl_puc_sumet = m_hlt_met_topocl_puc->sumEt() * 0.001;
  hlt_topocl_puc_phi = atan2(hlt_topocl_puc_ey, hlt_topocl_puc_ex);

  h_hlt_topocl_puc_ex->Fill( hlt_topocl_puc_ex ); // GeV
  h_hlt_topocl_puc_ey->Fill( hlt_topocl_puc_ey ); // GeV
  h_hlt_topocl_puc_met->Fill( hlt_topocl_puc_met ); // GeV
  h_hlt_topocl_puc_sumet->Fill( hlt_topocl_puc_sumet ); // GeV
  h_hlt_topocl_puc_phi->Fill( hlt_topocl_puc_phi ); // GeV



  //------------
  // MUONS
  //------------

  /// full copy 
  // get muon container of interest
  const xAOD::MuonContainer* m_muons(0);
  if ( !m_event->retrieve( m_muons, "Muons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// shallow copy for muon calibration and smearing tool
  // create a shallow copy of the muons container for MET building
  std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > muons_shallowCopy = xAOD::shallowCopyContainer( *m_muons );
  xAOD::MuonContainer* muonSC = muons_shallowCopy.first;

  // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
  // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
  // The method is defined in the header file xAODBase/IParticleHelpers.h
  bool setLinksMuon = xAOD::setOriginalObjectLink(*m_muons,*muonSC);
  if(!setLinksMuon) {
    Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
    return StatusCode::FAILURE;
  }

  ///////////////////
  // For MET study //
  ///////////////////

  // iterate over our shallow copy
  for (const auto& muon : *muonSC) { // C++11 shortcut
    //Info("execute()", "  original muon pt = %.2f GeV", ((*muonSC_itr)->pt() * 0.001)); // just to print out something
    //if(((*muonSC_itr)->eta()) > 2.5)
    //Info("execute()", "  muon eta = %.2f ", ((*muonSC_itr)->eta())); // just to print out something
    //Info("execute()", "  corrected muon pt = %.2f GeV", ((*muonSC_itr)->pt() * 0.001));
    passMuonSelection(*muon, eventInfo, primVertex);

  } // end for loop over shallow copied muons

  ///////////////////
  // For VBF study //
  ///////////////////

  /*
  // iterate over our shallow copy
  for (const auto& muon : *muonSC) { // C++11 shortcut
    // VBF Muon Selection
    passMuonVBF(*muon, eventInfo, primVertex);
    //Info("execute()", "  VBF muon pt = %.2f GeV", (muon->pt() * 0.001));
  } // end for loop over shallow copied muons
  */



  //------------
  // ELECTRONS
  //------------

  /// full copy 
  // get electron container of interest
  const xAOD::ElectronContainer* m_electrons(0);
  if ( !m_event->retrieve( m_electrons, "Electrons" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Electron container. Exiting." );
    return EL::StatusCode::FAILURE;
  }


  /// shallow copy for electron calibration tool
  // create a shallow copy of the electrons container for MET building
  std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* > elec_shallowCopy = xAOD::shallowCopyContainer( *m_electrons );
  xAOD::ElectronContainer* elecSC = elec_shallowCopy.first;

  // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
  // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
  // The method is defined in the header file xAODBase/IParticleHelpers.h
  bool setLinksElec = xAOD::setOriginalObjectLink(*m_electrons,*elecSC);
  if(!setLinksElec) {
    Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
    return StatusCode::FAILURE;
  }

  ///////////////////
  // For MET study //
  ///////////////////

  // iterate over our shallow copy
  for (const auto& electron : *elecSC) { // C++11 shortcut
    //Info("execute()", "  original electron pt = %.2f GeV", ((*elecSC_itr)->pt() * 0.001));
    //Info("execute()", "  corrected electron pt = %.2f GeV", ((*elecSC_itr)->pt() * 0.001));
    passElectronSelection(*electron, eventInfo, primVertex);

  } // end for loop over shallow copied electrons

  ///////////////////
  // For VBF study //
  ///////////////////

  /*
  // iterate over our shallow copy
  for (const auto& electron : *elecSC) { // C++11 shortcut
    // VBF Electron Selection
    passElectronVBF(*electron, eventInfo, primVertex);
    //Info("execute()", "  VBF electron pt = %.2f GeV", (electron->pt() * 0.001));
  } // end for loop over shallow copied electrons
  */





  //------------
  // PHOTONS
  //------------

  /// full copy 
  // get photon container of interest
  const xAOD::PhotonContainer* m_photons(0);
  if ( !m_event->retrieve( m_photons, "Photons" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Photon container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// shallow copy for photon calibration tool
  // create a shallow copy of the photons container for MET building
  std::pair< xAOD::PhotonContainer*, xAOD::ShallowAuxContainer* > phot_shallowCopy = xAOD::shallowCopyContainer( *m_photons );
  xAOD::PhotonContainer* photSC = phot_shallowCopy.first;

  // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
  // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
  // The method is defined in the header file xAODBase/IParticleHelpers.h
  bool setLinksPhoton = xAOD::setOriginalObjectLink(*m_photons,*photSC);
  if(!setLinksPhoton) {
    Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
    return StatusCode::FAILURE;
  }

  ///////////////////
  // For MET study //
  ///////////////////
  // iterate over our shallow copy
  for (const auto& photon : *photSC) { // C++11 shortcut
    //Info("execute()", "  original photon pt = %.2f GeV", ((*photSC_itr)->pt() * 0.001));
    //Info("execute()", "  corrected photon pt = %.2f GeV", ((*photSC_itr)->pt() * 0.001));
    passPhotonSelection(*photon, eventInfo);
  } // end for loop over shallow copied photons

  ///////////////////
  // For VBF study //
  ///////////////////

  /*
  // iterate over our shallow copy
  for (const auto& photon : *photSC) { // C++11 shortcut
    // VBF Tau Selection
    passPhotonVBF(*photon, eventInfo); 
  } // end for loop over shallow copied photons
  */



  //------------
  // TAUS
  //------------

  /// full copy 
  // get tau container of interest
  const xAOD::TauJetContainer* m_taus(0);
  if ( !m_event->retrieve( m_taus, "TauJets" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Tau container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  m_tauOverlappingElectronLLHDecorator->initializeEvent();

  /// shallow copy for tau calibration tool
  // create a shallow copy of the taus container for MET building
  std::pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* > tau_shallowCopy = xAOD::shallowCopyContainer( *m_taus );
  xAOD::TauJetContainer* tauSC = tau_shallowCopy.first;

  // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
  // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
  // The method is defined in the header file xAODBase/IParticleHelpers.h
  bool setLinksTau = xAOD::setOriginalObjectLink(*m_taus,*tauSC);
  if(!setLinksTau) {
    Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
    return StatusCode::FAILURE;
  }

  ///////////////////
  // For MET study //
  ///////////////////
  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    //Info("execute()", "  original tau pt = %.2f GeV", ((*tauSC_itr)->pt() * 0.001));
    //Info("execute()", "  corrected tau pt = %.2f GeV", ((*tauSC_itr)->pt() * 0.001));
    passTauSelection(*taujet, eventInfo);
  } // end for loop over shallow copied taus

  ///////////////////
  // For VBF study //
  ///////////////////

  /*
  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    // VBF Tau Selection
    passTauVBF(*taujet, eventInfo);
  } // end for loop over shallow copied taus
  */



  //------------
  // JETS
  //------------

  /// full copy 
  // get jet container of interest
  const xAOD::JetContainer* m_jets(0);
  if ( !m_event->retrieve( m_jets, jetType ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Jet container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// shallow copy for jet calibration tool
  // create a shallow copy of the jets container for MET building
  std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > jet_shallowCopy = xAOD::shallowCopyContainer( *m_jets );
  xAOD::JetContainer* jetSC = jet_shallowCopy.first;

  // Decorate objects with ElementLink to their originals -- this is needed to retrieve the contribution of each object to the MET terms.
  // You should make sure that you use the tag xAODBase-00-00-22, which is available from AnalysisBase-2.0.11.
  // The method is defined in the header file xAODBase/IParticleHelpers.h
  bool setLinksJet = xAOD::setOriginalObjectLink(*m_jets,*jetSC);
  if(!setLinksJet) {
    Error("execute()", "Failed to set original object links -- MET rebuilding cannot proceed.");
    return StatusCode::FAILURE;
  }

  // iterate over our shallow copy
  for (const auto& jets : *jetSC) { // C++11 shortcut
    //Info("execute()", "  original jet pt = %.2f GeV", jets->pt() * 0.001);

    // According to https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15

    // JES calibration
    if ( !m_jetCalibration->applyCalibration(*jets).isSuccess() ){
      Error("execute()", "Failed to apply calibration to Jet objects. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    float jetPt = jets->pt() * 0.001; /// GeV
    float jetEta = jets->eta();

    // JES correction
    if (!m_isData){
      if (jetPt > 15. ){
        if ( m_jetUncertaintiesTool->applyCorrection(*jets) != CP::CorrectionCode::Ok){ // apply correction and check return code
          Error("execute()", "Failed to apply JES correction to Jet objects. Exiting." );
          return EL::StatusCode::FAILURE;
        }
      }
    }

    // JER smearing
    if (!m_isData){
      if ( m_jerSmearingTool->applyCorrection(*jets) != CP::CorrectionCode::Ok){ // apply correction and check return code
        Error("execute()", "Failed to apply JER smearing. Exiting. ");
        return EL::StatusCode::FAILURE;
      }
    }

    // JVT Tool
    float newjvt = m_jvtagup->updateJvt(*jets);
    acc_jvt(*jets) = newjvt;

    //Info("execute()", "  corrected jet pt = %.2f GeV", jets->pt() * 0.001);
    //Info("execute()", "  updated jet jvt = %.2f ", newjvt);

    dec_baseline(*jets) = false;
    selectDec(*jets) = false; // To select objects for Overlap removal

    // pT cut
    if (jetPt > m_jetPtCut && jetEta < m_jetEtaCut) {
      dec_baseline(*jets) = true;
      selectDec(*jets) = true; // To select objects for Overlap removal
    }

  } // end for loop over shallow copied jets


  //----------------
  // Overlap Removal
  //----------------

  if (m_doORtool){
    //auto m_orTool = m_toolBox->getMasterHandle();
    if ( !m_orTool->removeOverlaps(elecSC, muonSC, jetSC, tauSC, photSC).isSuccess() ){
      Error("execute()", "Failed to apply the overlap removal to all objects. Exiting." );
      return EL::StatusCode::FAILURE;
    }


    // Now, dump all of the results

    // electrons
    for(auto electron : *elecSC){
      if(overlapAcc(*electron)) {
        nOverlapElectrons++;
        //Info("execute()", "  EventNumber : %i |  Overlap electron pt = %.2f GeV", EventNumber, (electron->pt() * 0.001));
      }
      nInputElectrons++;
    }
    // muons
    for(auto muon : *muonSC){
      if(overlapAcc(*muon)) nOverlapMuons++;
      nInputMuons++;
    }
    // jets
    for (auto jet : *jetSC) {
      if(overlapAcc(*jet)){
        nOverlapJets++;
        //Info("execute()", "  EventNumber : %i |  Overlap jet pt = %.2f GeV", EventNumber, (jet->pt() * 0.001));
      }
      nInputJets++;
    }
    // taus
    for(auto tau : *tauSC){
      if(overlapAcc(*tau)) nOverlapTaus++;
      nInputTaus++;
    }
    // photons
    for(auto photon : *photSC){
      if(overlapAcc(*photon)) nOverlapPhotons++;
      nInputPhotons++;
    }
  }


  // --------------------------------------------------------
  // Select Signal Jet and Bad Jet applying overlap removal
  // --------------------------------------------------------
  /// Creating New Hard Object Containers
  // [For jet identification] filter the Jet container m_jets, placing selected jets into m_goodJet
  xAOD::JetContainer* m_goodJet = new xAOD::JetContainer(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Jet>

  bool isBadJet = false;

  // iterate over our shallow copy
  for (const auto& jets : *jetSC) { // C++11 shortcut

    // Veto Jet (cleaning Jet)
    if (IsBadJet(*jets)) isBadJet = true;

    // Jet Signal Selection
    if (IsSignalJet(*jets)) {
      double jetPt = (jets->pt()) * 0.001; /// GeV
      h_jet_selection_pt->Fill( jetPt ); // GeV

      m_goodJet->push_back( jets );
    }
  } // end for loop over shallow copied jets



  //------------------------------------
  // Event Cleaning (Jet cleaning tool)
  //------------------------------------
  if (isBadJet) return EL::StatusCode::SUCCESS;
  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("Jet Cleaning");
  m_eventCutflow[4]+=1;



  //-----------------
  // Rebuild the MET
  //-----------------

  // Create a MissingETContainer with its aux store for each systematic
  xAOD::MissingETContainer* m_met = new xAOD::MissingETContainer();
  xAOD::MissingETAuxContainer* m_metAux = new xAOD::MissingETAuxContainer();
  m_met->setStore(m_metAux);

  //retrieve the original containers
  const xAOD::MissingETContainer* m_metCore(0);
  std::string coreMetKey = "MET_Core_" + jetType;
  coreMetKey.erase(coreMetKey.length() - 4); //this removes the Jets from the end of the jetType
  if ( !m_event->retrieve( m_metCore, coreMetKey ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Unable to retrieve MET core container: " );
    return EL::StatusCode::FAILURE;
  }

  //retrieve the MET association map
  const xAOD::MissingETAssociationMap* m_metMap(0);
  std::string metAssocKey = "METAssoc_" + jetType;
  metAssocKey.erase(metAssocKey.length() - 4 );//this removes the Jets from the end of the jetType
  if ( !m_event->retrieve( m_metMap, metAssocKey ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Unable to retrieve MissingETAssociationMap: " );
    return EL::StatusCode::FAILURE;
  }

  // do CST or TST
  //std::string softTerm = "SoftClus";
  std::string softTerm = "PVSoftTrk";

  // It is necessary to reset the selected objects before every MET calculation
  m_metMap->resetObjSelectionFlags();



  //==============================
  // For rebuild the emulated MET
  //==============================


  // Create a MissingETContainer with its aux store for each systematic
  xAOD::MissingETContainer* m_emulmet = new xAOD::MissingETContainer();
  xAOD::MissingETAuxContainer* m_emulmetAux = new xAOD::MissingETAuxContainer();
  m_emulmet->setStore(m_emulmetAux);

  //retrieve the original containers
  const xAOD::MissingETContainer* m_emulmetCore(0);
  if ( !m_event->retrieve( m_emulmetCore, coreMetKey ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Unable to retrieve MET core container: " );
    return EL::StatusCode::FAILURE;
  }

  //retrieve the MET association map
  const xAOD::MissingETAssociationMap* m_emulmetMap(0);
  if ( !m_event->retrieve( m_emulmetMap, metAssocKey ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Unable to retrieve MissingETAssociationMap: " );
    return EL::StatusCode::FAILURE;
  }

  // It is necessary to reset the selected objects before every MET calculation
  m_emulmetMap->resetObjSelectionFlags();





  //here we apply some basic cuts and rebuild the met at each step

  // Electron
  //-----------------

  /// Creat New Hard Object Containers
  // [For MET building] filter the Electron container m_electrons, placing selected electrons into m_MetElectrons
  ConstDataVector<xAOD::ElectronContainer> m_MetElectrons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Electron>

  // iterate over our shallow copy
  for (const auto& electron : *elecSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*electron)) {
      if(m_doORtool){
        if(!overlapAcc(*electron)) {
          //double elecPt = (electron->pt()) * 0.001; /// GeV
          //Info("execute()", "  Selected MET electron pt from shallow copy = %.2f GeV", ((*elecSC_itr)->pt() * 0.001));
          //Info("execute()", "  Selected MET electron pt from new Electron Container = %.2f GeV", (electron->pt() * 0.001));  
          m_MetElectrons.push_back( electron );
        }
      }
      else m_MetElectrons.push_back( electron );
    }
  } // end for loop over shallow copied electrons
  //const xAOD::ElectronContainer* p_MetElectrons = m_MetElectrons.asDataVector();
  m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
      xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
      m_met,                                      //filling this met container
      m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
      m_metMap);                                  //and this association map

  if (m_isZee || m_isWenu){
    m_metMaker->markInvisible(m_MetElectrons.asDataVector(), m_emulmetMap);
  }
  else if (m_isZmumu || m_isWmunu){
    m_metMaker->rebuildMET("RefElectron",           //name of metElectrons in metContainer
        xAOD::Type::Electron,                       //telling the rebuilder that this is electron met
        m_emulmet,                                  //filling this met container
        m_MetElectrons.asDataVector(),              //using these metElectrons that accepted our cuts
        m_emulmetMap);                              //and this association map
  }



  // Photon
  //-----------------

  /// Creat New Hard Object Containers
  // [For MET building] filter the Photon container m_photons, placing selected photons into m_MetPhotons
  ConstDataVector<xAOD::PhotonContainer> m_MetPhotons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Photon>

  // iterate over our shallow copy
  for (const auto& photon : *photSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*photon)) {
      if(m_doORtool){
        if(!overlapAcc(*photon)){
          //double photPt = (photon->pt()) * 0.001; /// GeV
          //Info("execute()", "  Selected MET photon pt from shallow copy = %.2f GeV", ((*photSC_itr)->pt() * 0.001));
          //Info("execute()", "  Selected MET photon pt from new Photon Container = %.2f GeV", (phot->pt() * 0.001));  
          m_MetPhotons.push_back( photon );
        }
      }
      else m_MetPhotons.push_back( photon );
    }
  } // end for loop over shallow copied photons
  m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
      xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
      m_met,                                    //filling this met container
      m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
      m_metMap);                                //and this association map

  if (m_isZmumu || m_isZee || m_isWmunu || m_isWenu){
    m_metMaker->rebuildMET("RefPhoton",           //name of metPhotons in metContainer
        xAOD::Type::Photon,                       //telling the rebuilder that this is photon met
        m_emulmet,                                //filling this met container
        m_MetPhotons.asDataVector(),              //using these metPhotons that accepted our cuts
        m_emulmetMap);                            //and this association map
  }

  // TAUS
  //-----------------

  /// Creat New Hard Object Containers
  // [For MET building] filter the TauJet container m_taus, placing selected taus into m_MetTaus
  ConstDataVector<xAOD::TauJetContainer> m_MetTaus(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::TauJet>

  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*taujet)) {
      if(m_doORtool){
        if(!overlapAcc(*taujet)){
          //double tauPt = (taujet->pt()) * 0.001; /// GeV
          //Info("execute()", "  Selected MET tau pt from shallow copy = %.2f GeV", ((*tauSC_itr)->pt() * 0.001));
          //Info("execute()", "  Selected MET tau pt from new Tau Container = %.2f GeV", (tau->pt() * 0.001));  
          m_MetTaus.push_back( taujet );
        }
      }
      else m_MetTaus.push_back( taujet );
    }
  } // end for loop over shallow copied taus
  m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
      xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
      m_met,                                 //filling this met container
      m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
      m_metMap);                             //and this association map

  if (m_isZmumu || m_isZee || m_isWmunu || m_isWenu){
    m_metMaker->rebuildMET("RefTau",           //name of metTaus in metContainer
        xAOD::Type::Tau,                       //telling the rebuilder that this is tau met
        m_emulmet,                             //filling this met container
        m_MetTaus.asDataVector(),              //using these metTaus that accepted our cuts
        m_emulmetMap);                         //and this association map
  }

  // Muon
  //-----------------

  /// Creat New Hard Object Containers
  // [For MET building] filter the Muon container m_muons, placing selected muons into m_MetMuons
  ConstDataVector<xAOD::MuonContainer> m_MetMuons(SG::VIEW_ELEMENTS); // This is really a DataVector<xAOD::Muon>

  // iterate over our shallow copy
  for (const auto& muon : *muonSC) { // C++11 shortcut
    // For MET rebuilding
    if (dec_baseline(*muon)) {
      if(m_doORtool){
        if(!overlapAcc(*muon)){
          //double muPt = (muon->pt()) * 0.001; /// GeV
          //Info("execute()", "  Selected Muon pt from shallow copy = %.2f GeV", ((*muonSC_itr)->pt() * 0.001));
          //Info("execute()", "  Selected muon pt from new Muon Container = %.2f GeV", (muon->pt() * 0.001));  
          m_MetMuons.push_back( muon );
        }
      }
      else m_MetMuons.push_back( muon );
    }
  } // end for loop over shallow copied muons
  m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
      xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
      m_met,                                  //filling this met container
      m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
      m_metMap);                              //and this association map

  if (m_isZmumu || m_isWmunu){
    m_metMaker->markInvisible(m_MetMuons.asDataVector(), m_emulmetMap);
    met::addGhostMuonsToJets(*m_muons, *jetSC);
  }
  else if (m_isZee || m_isWenu){
    m_metMaker->rebuildMET("RefMuon",           //name of metMuons in metContainer
        xAOD::Type::Muon,                       //telling the rebuilder that this is muon met
        m_emulmet,                              //filling this met container
        m_MetMuons.asDataVector(),              //using these metMuons that accepted our cuts
        m_emulmetMap);                          //and this association map
  }

  // JET
  //-----------------

  //Now time to rebuild jetMet and get the soft term
  //This adds the necessary soft term for both CST and TST
  //these functions create an xAODMissingET object with the given names inside the container
  m_metMaker->rebuildJetMET("RefJet",          //name of jet met
      "SoftClus",        //name of soft cluster term met
      "PVSoftTrk",       //name of soft track term met
      m_met,             //adding to this new met container
      jetSC,             //using this jet collection to calculate jet met
      m_metCore,         //core met container
      m_metMap,          //with this association map
      true);             //apply jet jvt cut

  if (m_isZmumu || m_isZee || m_isWmunu || m_isWenu){
    m_metMaker->rebuildJetMET("RefJet",          //name of jet met
        "SoftClus",        //name of soft cluster term met
        "PVSoftTrk",       //name of soft track term met
        m_emulmet,         //adding to this new met container
        jetSC,             //using this jet collection to calculate jet met
        m_emulmetCore,     //core met container
        m_emulmetMap,      //with this association map
        true);             //apply jet jvt cut
  }

  // MET Build
  //-----------------

  //m_metMaker->rebuildTrackMET("RefJetTrk", softTerm, m_met, jetSC, m_metCore, m_metMap, true);

  //this builds the final track or cluster met sums, using systematic varied container
  //In the future, you will be able to run both of these on the same container to easily output CST and TST
  m_metMaker->buildMETSum("Final", m_met, (*m_met)[softTerm]->source());

  if (m_isZmumu || m_isZee || m_isWmunu || m_isWenu){
    m_metMaker->buildMETSum("Final", m_emulmet, (*m_emulmet)[softTerm]->source());
  }


  // Fill MET RefFinal
  float MET_ex = -9e9;
  float MET_ey = -9e9;
  float MET = -9e9;
  float SumET = -9e9;
  float MET_phi = -9e9;

  MET_ex = ((*m_met)["Final"]->mpx()) * 0.001;
  MET_ey = ((*m_met)["Final"]->mpy()) * 0.001;
  //MET = sqrt(MET_ex*MET_ex+MET_ey*MET_ey);
  MET = ((*m_met)["Final"]->met()) * 0.001;
  SumET = ((*m_met)["Final"]->sumet()) * 0.001;
  //MET_phi = atan2(MET_ey, MET_ex);
  MET_phi = ((*m_met)["Final"]->phi());

  h_met_ex->Fill( MET_ex ); // GeV
  h_met_ey->Fill( MET_ey ); // GeV
  h_met->Fill( MET ); // GeV
  h_sumet->Fill( SumET ); // GeV
  h_met_phi->Fill( MET_phi ); // GeV


  // Fill emulated MET RefFinal
  float emulMET_ex = -9e9;
  float emulMET_ey = -9e9;
  float emulMET = -9e9;
  float emulSumET = -9e9;
  float emulMET_phi = -9e9;

  if (m_isZmumu || m_isZee || m_isWmunu || m_isWenu){
    emulMET_ex = ((*m_emulmet)["Final"]->mpx()) * 0.001;
    emulMET_ey = ((*m_emulmet)["Final"]->mpy()) * 0.001;
    emulMET = ((*m_emulmet)["Final"]->met()) * 0.001;
    emulSumET = ((*m_emulmet)["Final"]->sumet()) * 0.001;
    emulMET_phi = ((*m_emulmet)["Final"]->phi());

    h_emulmet_ex->Fill( emulMET_ex ); // GeV
    h_emulmet_ey->Fill( emulMET_ey ); // GeV
    h_emulmet->Fill( emulMET ); // GeV
    h_emulsumet->Fill( emulSumET ); // GeV
    h_emulmet_phi->Fill( emulMET_phi ); // GeV
  }



  //-----------------------------------------------
  // Define Good Leptons and Calculate Scale Factor
  //-----------------------------------------------

  ///////////////
  // Good Muon //
  ///////////////
  xAOD::MuonContainer* m_goodMuon = new xAOD::MuonContainer(SG::VIEW_ELEMENTS);
  // iterate over our shallow copy
  for (const auto& muon : *muonSC) { // C++11 shortcut
    // Muon Selection
    //if (dec_baseline(*muon)) { // For VBF study
    if (dec_signal(*muon)) { // For MET Trigger study
      if(m_doORtool){
        if (!overlapAcc(*muon)){
          m_goodMuon->push_back( muon );
        }
      }
      else m_goodMuon->push_back( muon );
      //Info("execute()", "  Good muon pt = %.2f GeV", (muon->pt() * 0.001));
    }
  } // end for loop over shallow copied muons

  ///////////////////
  // Good Electron //
  ///////////////////
  xAOD::ElectronContainer* m_goodElectron = new xAOD::ElectronContainer(SG::VIEW_ELEMENTS);
  // iterate over our shallow copy
  for (const auto& electron : *elecSC) { // C++11 shortcut
    // Electron Selection
    //if (dec_baseline(*electron)) { // For VBF study
    if (dec_signal(*electron)) { // MET Trigger sutdy
      if(m_doORtool){
        if(!overlapAcc(*electron)){
          m_goodElectron->push_back( electron );
        }
      }
      else m_goodElectron->push_back( electron );
      //Info("execute()", "  Good electron pt = %.2f GeV", (electron->pt() * 0.001));
    }
  } // end for loop over shallow copied electrons

  //////////////
  // Good Tau //
  //////////////
  xAOD::TauJetContainer* m_goodTau = new xAOD::TauJetContainer(SG::VIEW_ELEMENTS);
  // iterate over our shallow copy
  for (const auto& taujet : *tauSC) { // C++11 shortcut
    // Tau Selection
    if (dec_baseline(*taujet)) {
      if(m_doORtool){
        if(!overlapAcc(*taujet)){
          m_goodTau->push_back( taujet );
        }
      }
      else m_goodTau->push_back( taujet );
      //Info("execute()", "  Good tau pt = %.2f GeV", (taujet->pt() * 0.001));
    }
  } // end for loop over shallow copied electrons


  /////////////////////////
  // Overlap removed Jet //
  /////////////////////////
  xAOD::JetContainer* m_signalJet = new xAOD::JetContainer(SG::VIEW_ELEMENTS);
  // iterate over our shallow copy
  for (const auto& jet : *m_goodJet) { // C++11 shortcut

    if (m_doORmanual) {
      bool isORjet = false;

      for (const auto& muon : *m_goodMuon) {
        if (DeltaR(jet->eta(), muon->eta(), jet->phi(), muon->phi()) < m_ORJETdeltaR) isORjet = true;
      }

      for (const auto& electron : *m_goodElectron) {
        if (DeltaR(jet->eta(), electron->eta(), jet->phi(), electron->phi()) < m_ORJETdeltaR) isORjet = true;
      }

      for (const auto& tau : *m_goodTau) {
        if (DeltaR(jet->eta(), tau->eta(), jet->phi(), tau->phi()) < m_ORJETdeltaR) isORjet = true;
      }

      if (!isORjet) m_signalJet->push_back( jet );
    }
    else m_signalJet->push_back( jet );

  } // end for loop over shallow copied jets


  //------------------------
  // Define DiJet Properties
  // -----------------------

  TLorentzVector jet1;
  TLorentzVector jet2;
  float jet1_pt = 0;
  float jet2_pt = 0;
  float jet1_phi = 0;
  float jet2_phi = 0;
  float jet1_rapidity = 0;
  float jet2_rapidity = 0;

  float signalJet_ht = 0;
  float dPhiJet1Met = 0;
  float dPhiJet2Met = 0;

  float mjj = 0;
  bool pass_diJet = false; // Select DiJet
  bool pass_CJV = true; // Central Jet Veto (CJV)
  bool pass_dPhijetmet = true; // deltaPhi(Jet,MET)
  bool pass_dPhiDijetMet = true; // deltaPhi(Jet1,MET) or deltaPhi(Jet2,MET)


  // DiJet Selection
  if (m_signalJet->size() > 1) {
    std::partial_sort(m_signalJet->begin(), m_signalJet->begin()+2, m_signalJet->end(), DescendingPt());

    jet1 = m_signalJet->at(0)->p4();
    jet2 = m_signalJet->at(1)->p4();
    jet1_pt = m_signalJet->at(0)->pt() * 0.001;
    jet2_pt = m_signalJet->at(1)->pt() * 0.001;
    jet1_phi = m_signalJet->at(0)->phi();
    jet2_phi = m_signalJet->at(1)->phi();
    jet1_rapidity = m_signalJet->at(0)->rapidity();
    jet2_rapidity = m_signalJet->at(1)->rapidity();
    auto dijet = jet1 + jet2;
    mjj = dijet.M() * 0.001;

    //Info("execute()", "  jet1 = %.2f GeV, jet2 = %.2f GeV", jet1_pt, jet2_pt);
    //Info("execute()", "  mjj = %.2f GeV", mjj);

    if ( jet1_pt >  m_diJet1PtCut && jet2_pt > m_diJet2PtCut ){
      if ( jet1_rapidity < m_diJetEtaCut && jet2_rapidity < m_diJetEtaCut ){
        if ( m_jetCleaningTight->accept( *m_signalJet->at(0) ) ){ //Tight Leading Jet 
          pass_diJet = true;
        }
      }
    }

    // deltaPhi(Jet1,MET) or deltaPhi(Jet2,MET) decision
    if (m_isZmumu || m_isZee || m_isWmunu || m_isWenu){
      dPhiJet1Met = DeltaPhi(jet1_phi, emulMET_phi);
      dPhiJet2Met = DeltaPhi(jet2_phi, emulMET_phi);
    }
    else {
      dPhiJet1Met = DeltaPhi(jet1_phi, MET_phi);
      dPhiJet2Met = DeltaPhi(jet2_phi, MET_phi);
    }
    if ( dPhiJet1Met < 0.4 ) pass_dPhiDijetMet = false ;
    if ( dPhiJet2Met < 0.4 ) pass_dPhiDijetMet = false ;

  } // DiJet selection loop


  // N_jet >= 1
  if (m_signalJet->size() > 0) {

    // loop over the jets in the Good Jets Container
    for (const auto& jet : *m_signalJet) {
      float signal_jet_pt = jet->pt() * 0.001;
      float signal_jet_rapidity = jet->rapidity();
      float signal_jet_phi = jet->phi();

      // dphijetmet
      float dPhijetmet = DeltaPhi(signal_jet_phi,MET_phi);
      //Info("execute()", "  delta phi = %.2f", dPhijetmet);
      if ( dPhijetmet < 0.5 ) pass_dPhijetmet = false;

      // Central Jet Veto (CJV)
      if ( m_signalJet->size() > 2 && pass_diJet ){
        if (m_signalJet->at(0) != jet && m_signalJet->at(1) != jet){
          //cout << "m_signalJet->at(0) = " << m_signalJet->at(0) << " jet = " << jet << endl;
          if (signal_jet_pt > m_CJVptCut) {
            if ( (jet1_rapidity > jet2_rapidity) && (signal_jet_rapidity < jet1_rapidity && signal_jet_rapidity > jet2_rapidity)){
              pass_CJV = false;
            }
            if ( (jet1_rapidity < jet2_rapidity) && (signal_jet_rapidity > jet1_rapidity && signal_jet_rapidity < jet2_rapidity)){
              pass_CJV = false;
            }
          }
        }
      }

      //Info("execute()", "  Zvv Signal Jet pt = %.2f GeV, eta = %.2f", signal_pt_jet, signal_eta_jet);
      signalJet_ht += signal_jet_pt;
    } // Jet loop

  } //N_jet >= 1

 

  //----------------------------------
  // Define Zmumu and Wmunu Selection
  //----------------------------------

  bool pass_Zmumu = false; // Select Zmumu channel
  bool pass_Wmunu = false; // Select Wmunu channel
  float mll_muon = 0.; // For Zmumu channel
  float mT_muon = 0.;// For Wmunu channel
  if (m_isZmumu || m_isWmunu){
    // Zmumu muons
    TLorentzVector muon1;
    TLorentzVector muon2;
    float muon1_pt = 0.;
    float muon2_pt = 0.;
    float muon1_charge = 0.;
    float muon2_charge = 0.;
    // Wmunu muon
    float muon_pt = 0.;
    float muon_phi = 0.;

    // Zmumu Selection
    if (m_goodMuon->size() > 1) {
      std::partial_sort(m_goodMuon->begin(), m_goodMuon->begin()+2, m_goodMuon->end(), DescendingPt());

      muon1 = m_goodMuon->at(0)->p4();
      muon2 = m_goodMuon->at(1)->p4();
      muon1_pt = m_goodMuon->at(0)->pt() * 0.001;
      muon2_pt = m_goodMuon->at(1)->pt() * 0.001;
      muon1_charge = m_goodMuon->at(0)->charge();
      muon2_charge = m_goodMuon->at(1)->charge();
      auto Zmass_muon = muon1 + muon2;
      mll_muon = Zmass_muon.M() * 0.001;

      //Info("execute()", "  muon1 = %.2f GeV, muon2 = %.2f GeV", muon1_pt, muon2_pt);
      //Info("execute()", "  mll (Zmumu) = %.2f GeV", mll_muon);

      if ( muon1_pt >  m_LeadLepPtCut && muon2_pt > m_SubLeadLepPtCut ){
        if ( muon1_charge * muon2_charge < 0 ){
          pass_Zmumu = true;
          //Info("execute()", "  Leading muon = %.2f GeV, Subleading muon = %.2f GeV", muon1_pt, muon2_pt);
        }
      }
    } // Zmumu selection loop

    // Wmunu Selection
    if (m_goodMuon->size() == 1) {
      muon_pt = m_goodMuon->at(0)->pt() * 0.001;
      muon_phi = m_goodMuon->at(0)->phi();
      mT_muon = TMath::Sqrt( 2. * muon_pt * MET * ( 1. - TMath::Cos(muon_phi - MET_phi) ) );

      if ( muon_pt > 25. ){
        pass_Wmunu = true;
      }
    } //Wmunu Selection loop

  }


  //----------------------
  // Define Zee Selection
  // ---------------------

  bool pass_Zee = false; // Select Zee channel
  float mll_electron = 0;
  if (m_isZee){
    TLorentzVector electron1;
    TLorentzVector electron2;
    float electron1_pt = 0;
    float electron2_pt = 0;
    float electron1_charge = 0;
    float electron2_charge = 0;



    // Zee Selection
    if (m_goodElectron->size() > 1) {
      std::partial_sort(m_goodElectron->begin(), m_goodElectron->begin()+2, m_goodElectron->end(), DescendingPt());

      electron1 = m_goodElectron->at(0)->p4();
      electron2 = m_goodElectron->at(1)->p4();
      electron1_pt = m_goodElectron->at(0)->pt() * 0.001;
      electron2_pt = m_goodElectron->at(1)->pt() * 0.001;
      electron1_charge = m_goodElectron->at(0)->charge();
      electron2_charge = m_goodElectron->at(1)->charge();
      auto Zmass_electron = electron1 + electron2;
      mll_electron = Zmass_electron.M() * 0.001;

      //Info("execute()", "  electron1 = %.2f GeV, electron2 = %.2f GeV", electron1_pt, electron2_pt);
      //Info("execute()", "  mll (Zee) = %.2f GeV", mll_electron);

      if ( electron1_pt >  m_LeadLepPtCut && electron2_pt > m_SubLeadLepPtCut ){
        if ( electron1_charge * electron2_charge < 0 ){
          pass_Zee = true;
          //Info("execute()", "  Leading electron = %.2f GeV, Subleading electron = %.2f GeV", electron1_pt, electron2_pt);
        }
      }

    } // Zee selection loop
  }




  // ------------------
  // Get isolated track
  // ------------------
/*
  // Retrieve main TrackParticle collection
  const xAOD::TrackParticleContainer* inTracks(0);
  if ( !m_event->retrieve( inTracks, "InDetTrackParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve TrackParticle container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  int Nisotrk = NumIsoTracks(inTracks, primVertex, 3., 10.) - NumMuonIsoTrack(muonSC, inTracks, primVertex, 3., 10.) - NumElecIsoTrack(elecSC, inTracks, primVertex, 3., 10.);

  bool passIsoTrk = true;
  if (Nisotrk > 0) {
    passIsoTrk = false;
    //Info("execute()", "  The number of Isolated track counted = %i (N_SignalMuon = %lu, N_SignalElec = %lu)", Nisotrk, m_signalMuon->size(), m_signalElectron->size() );
  }
*/



  //-----------
  // VBF study 
  //-----------

  //-------------------------------
  // Z -> vv + JET EVENT SELECTION
  //-------------------------------

  if (m_isZvv){
    if ( m_trigDecisionTool->isPassed("HLT_xe80") ) {
      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]MET Trigger");
      m_eventCutflow[5]+=1;
      if ( MET > m_metCut ) {
        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]MET cut");
        m_eventCutflow[6]+=1;
/*
        Info("execute()", "=====================================");
        Info("execute()", " Event # = %llu", eventInfo->eventNumber());
        Info("execute()", " Good Event number = %i", m_eventCutflow[6]);
        Info("execute()", " MET = %.3f GeV", MET);
        Info("execute()", " RefElectron = %.3f GeV", ((*m_met)["RefElectron"]->met()) * 0.001);
        Info("execute()", " RefPhoton = %.3f GeV", ((*m_met)["RefPhoton"]->met()) * 0.001);
        Info("execute()", " RefTau = %.3f GeV", ((*m_met)["RefTau"]->met()) * 0.001);
        Info("execute()", " RefMuon = %.3f GeV", ((*m_met)["RefMuon"]->met()) * 0.001);
        Info("execute()", " RefJet = %.3f GeV", ((*m_met)["RefJet"]->met()) * 0.001);
        Info("execute()", " SoftClus = %.3f GeV", ((*m_met)["SoftClus"]->met()) * 0.001);
        Info("execute()", " PVSoftTrk = %.3f GeV", ((*m_met)["PVSoftTrk"]->met()) * 0.001);
        Info("execute()", " # of good jets = %lu", m_signalJet->size());
        if (m_signalJet->size() > 0){
          int jetCount = 0;
          for (const auto& jet : *m_signalJet) {
            jetCount++;
            Info("execute()", " jet # : %i", jetCount);
            Info("execute()", " jet pt = %.3f GeV", jet->pt() * 0.001);
            Info("execute()", " jet eta = %.3f GeV", jet->eta());
            Info("execute()", " jet phi = %.3f GeV", jet->phi());
          }
        }
        if (m_goodElectron->size() > 0){
          int eleCount = 0;
          for (const auto& electron : *m_goodElectron) {
            eleCount++;
            Info("execute()", " electron # : %i", eleCount);
            Info("execute()", " electron pt = %.3f GeV", electron->pt() * 0.001);
            Info("execute()", " electron eta = %.3f GeV", electron->eta());
            Info("execute()", " electron phi = %.3f GeV", electron->phi());
          }
        }
        if (m_goodMuon->size() > 0){
          int muCount = 0;
          for (const auto& muon : *m_goodMuon) {
            muCount++;
            Info("execute()", " muon # : %i", muCount);
            Info("execute()", " muon pt = %.3f GeV", muon->pt() * 0.001);
            Info("execute()", " muon eta = %.3f GeV", muon->eta());
            Info("execute()", " muon phi = %.3f GeV", muon->phi());
          }
        }
        if (m_goodTau->size() > 0){
          int tauCount = 0;
          for (const auto& tau : *m_goodTau) {
            tauCount++;
            Info("execute()", " tau # : %i", tauCount);
            Info("execute()", " tau pt = %.3f GeV", tau->pt() * 0.001);
            Info("execute()", " tau eta = %.3f GeV", tau->eta());
            Info("execute()", " tau phi = %.3f GeV", tau->phi());
          }
        }
*/
        if (m_goodElectron->size() == 0) {
          if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]Electron Veto");
          m_eventCutflow[7]+=1;
          if ( m_goodMuon->size() == 0) {
            if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]Muon Veto");
            m_eventCutflow[8]+=1;
            if (m_goodTau->size() == 0) {
              if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]Tau Veto");
              m_eventCutflow[9]+=1;
              if ( m_signalJet->size() > 1 ) {
                if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]At least two Jets");
                m_eventCutflow[10]+=1;
                if ( pass_diJet ) {
                  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]DiJet");
                  m_eventCutflow[11]+=1;
                  if ( mjj > m_mjjCut ) {
                    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]mjj cut");
                    m_eventCutflow[12]+=1;
                    if ( pass_CJV ) {
                      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]CJV");
                      m_eventCutflow[13]+=1;
                      if ( pass_dPhiDijetMet ) {
                        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zvv]dPhi(j,MET) cut");
                        m_eventCutflow[14]+=1;
                        //h_zvv_offline_met->Fill( MET ); // GeV
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  //---------------------------------
  // Z -> mumu + JET EVENT SELECTION
  //---------------------------------

  if (m_isZmumu){
    if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("MET Trigger");
      m_eventCutflow[5]+=1;
      if ( emulMET > m_metCut ) {
        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]MET cut");
        m_eventCutflow[6]+=1;
        if (m_goodElectron->size() == 0) {
          if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]Electron Veto");
          m_eventCutflow[7]+=1;
          if ( m_goodMuon->size() > 1) {
            if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]At least Two Electrons");
            m_eventCutflow[8]+=1;
            if (m_goodTau->size() == 0) {
              if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]Tau Veto");
              m_eventCutflow[9]+=1;
              if ( pass_Zmumu && m_goodMuon->size() == 2 && mll_muon > 66. && mll_muon < 116. ){
                if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]mll cut");
                m_eventCutflow[10]+=1;
                if ( m_signalJet->size() > 1 ) {
                  if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]At least Two Jets");
                  m_eventCutflow[11]+=1;
                  if ( pass_diJet ) {
                    if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]DiJet");
                    m_eventCutflow[12]+=1;
                    if ( mjj > m_mjjCut ) {
                      if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]mjj cut");
                      m_eventCutflow[13]+=1;
                      if ( pass_CJV ) {
                        if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]CJV cut");
                        m_eventCutflow[14]+=1;
                        if ( pass_dPhiDijetMet ) {
                          if (m_useBitsetCutflow) m_BitsetCutflow->FillCutflow("[Zmumu]dPhi(j,MET) cut");
                          m_eventCutflow[15]+=1;
                          // Calculate muon SF
                          if (!m_isData) {
                            double totalMuonSF_Zmumu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoSF, m_ttvaSF);
                            //Info("execute()", " Zmumu Total Muon SF = %.3f ", totalMuonSF_Zmumu);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  //---------------------------------
  // W -> munu + JET EVENT SELECTION
  //---------------------------------

  if (m_isWmunu){
    if ( m_trigDecisionTool->isPassed("HLT_xe70") ) {
      if ( emulMET > m_metCut ) {
        if (m_goodElectron->size() == 0) {
          if ( m_goodMuon->size() > 0 ) {
            if (m_goodTau->size() == 0) {
              if ( pass_Wmunu && m_goodMuon->size() == 1 && mT_muon > 30. && mT_muon < 100. ){
                if ( m_signalJet->size() > 1 ) {
                  if ( pass_diJet ) {
                    if ( mjj > m_mjjCut ) {
                      if ( pass_CJV ) {
                        if ( pass_dPhiDijetMet ) {
                          // Calculate muon SF
                          if (!m_isData) {
                            double totalMuonSF_Wmunu = GetTotalMuonSF(*m_goodMuon, m_recoSF, m_isoSF, m_ttvaSF);
                            Info("execute()", " Wmunu Total Muon SF = %.3f ", totalMuonSF_Wmunu);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }




  //-------------------------------
  // Z -> ee + JET EVENT SELECTION
  //-------------------------------

  if (m_isZee){
    if ( m_trigDecisionTool->isPassed("HLT_e24_lhmedium_L1EM20VH") || m_trigDecisionTool->isPassed("HLT_e60_lhmedium") || m_trigDecisionTool->isPassed("HLT_e120_lhloose") ) {
    //if (m_trigDecisionTool->isPassed("HLT_e24_lhmedium_iloose_L1EM20VH") || m_trigDecisionTool->isPassed("HLT_e60_lhmedium*")) {
      m_eventCutflow[5]+=1;
      if ( emulMET > m_metCut ) {
        m_eventCutflow[6]+=1;
        if (m_goodElectron->size() > 1) {
          m_eventCutflow[7]+=1;
          if ( m_goodMuon->size() == 0) {
            m_eventCutflow[8]+=1;
            if (m_goodTau->size() == 0) {
              m_eventCutflow[9]+=1;
              if ( pass_Zee && m_goodElectron->size() == 2 && mll_electron > 66. && mll_electron < 116. ) {
                m_eventCutflow[10]+=1;
                if ( m_signalJet->size() > 1 ) {
                  m_eventCutflow[11]+=1;
                  if ( pass_diJet ) {
                    m_eventCutflow[12]+=1;
                    if ( mjj > m_mjjCut ) {
                      m_eventCutflow[13]+=1;
                      if ( pass_CJV ) {
                        m_eventCutflow[14]+=1;
                        if ( pass_dPhiDijetMet ) {
                          m_eventCutflow[15]+=1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }





  //////////////////////////
  // Delete copy containers
  //////////////////////////

  // Deep copies. Clearing containers deletes contents including AuxStore.
  delete m_goodJet;

  // MET study
  //delete m_signalMuon;
  //delete m_signalElectron;

  // VBF study
  delete m_goodMuon;
  delete m_goodElectron;
  delete m_goodTau;
  delete m_signalJet;


  //////////////////////////////////
  // Delete shallow copy containers
  //////////////////////////////////

  // The containers created by the shallow copy are owned by you. Remember to delete them
  delete m_met;
  delete m_metAux;

  delete m_emulmet;
  delete m_emulmetAux;

  delete muons_shallowCopy.first;
  delete muons_shallowCopy.second;

  delete elec_shallowCopy.first;
  delete elec_shallowCopy.second;

  delete phot_shallowCopy.first;
  delete phot_shallowCopy.second;

  delete tau_shallowCopy.first;
  delete tau_shallowCopy.second;

  delete jet_shallowCopy.first;
  delete jet_shallowCopy.second;



  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.


  // cutflow
  if (m_useBitsetCutflow) m_BitsetCutflow->PushBitSet();;

  //*************************
  // deleting of all tools
  // ************************

  // GRL
  if (m_grl) {
    delete m_grl;
    m_grl = 0;
  }

  // cleaning up trigger tools
  if( m_trigConfigTool ) {
    delete m_trigConfigTool;
    m_trigConfigTool = 0;
  }
  if( m_trigDecisionTool ) {
    delete m_trigDecisionTool;
    m_trigDecisionTool = 0;
  }


  // in finalize, delete tool:
  /// Muon Calibration and Smearing Tool
  if(m_muonCalibrationAndSmearingTool){
    delete m_muonCalibrationAndSmearingTool;
    m_muonCalibrationAndSmearingTool = 0;
  }

  /// Electron Calibration and Smearing Tool
  if(m_egammaCalibrationAndSmearingTool){
    delete m_egammaCalibrationAndSmearingTool;
    m_egammaCalibrationAndSmearingTool = 0;
  }

  /// Muon selector tool
  if(m_muonSelection){
    delete m_muonSelection;
    m_muonSelection = 0;
  }

  /// Loose Muon selector tool
  if(m_loosemuonSelection){
    delete m_loosemuonSelection;
    m_loosemuonSelection = 0;
  }

  /// Electron selector tool
  if(m_LHToolTight2015){
    delete m_LHToolTight2015;
    m_LHToolTight2015 = 0;
  }
  if(m_LHToolMedium2015){
    delete m_LHToolMedium2015;
    m_LHToolMedium2015 = 0;
  }
  if(m_LHToolLoose2015){
    delete m_LHToolLoose2015;
    m_LHToolLoose2015 = 0;
  }


  /// Recomputing the photon ID flags
  if(m_photonTightIsEMSelector){
    delete m_photonTightIsEMSelector;
    m_photonTightIsEMSelector = 0;
  }
  if(m_photonMediumIsEMSelector){
    delete m_photonMediumIsEMSelector;
    m_photonMediumIsEMSelector = 0;
  }
  if(m_photonLooseIsEMSelector){
    delete m_photonLooseIsEMSelector;
    m_photonLooseIsEMSelector = 0;
  }

  /// IsolationSelectionTool
  if(m_IsolationSelectionTool){
    delete m_IsolationSelectionTool;
    m_IsolationSelectionTool = 0;
  }

  /// IsolationSelectionTool for VBF signal
  if(m_IsoToolVBF){
    delete m_IsoToolVBF;
    m_IsoToolVBF = 0;
  }

  /// Tau Smearing Tool
  if(m_tauSmearingTool){
    delete m_tauSmearingTool;
    m_tauSmearingTool = 0;
  }

  /// Tau Selection Tool
  if(m_tauSelTool){
    delete m_tauSelTool;
    m_tauSelTool = 0;
  }

  /// Tau Selection Tool for VBF signal
  if(m_tauSelToolVBF){
    delete m_tauSelToolVBF;
    m_tauSelToolVBF = 0;
  }

  /// JES Calibration
  if(m_jetCalibration){
    delete m_jetCalibration;
    m_jetCalibration = 0;
  }

  /// JES uncertainty
  if(m_jetUncertaintiesTool){
    delete m_jetUncertaintiesTool;
    m_jetUncertaintiesTool = 0;
  }

  /// JER Tool
  if(m_jerTool){
    delete m_jerTool;
    m_jerTool = 0;
  }

  ///  JER Smearing Tool
  if(m_jerSmearingTool){
    delete m_jerSmearingTool;
    m_jerSmearingTool = 0;
  }

  /// JVT Tool
  if(m_jvtag){
    delete m_jvtag;
    m_jvtag = 0;
  }

  /// Jet Cleaning Tool
  if(m_jetCleaningTight) {
    delete m_jetCleaningTight;
    m_jetCleaningTight = 0;
  }
  if(m_jetCleaningLoose) {
    delete m_jetCleaningLoose;
    m_jetCleaningLoose = 0;
  }

  /// MET Tool
  if(m_metMaker){
    delete m_metMaker;
    m_metMaker = 0;
  }

  /// Muon Efficiency Tool
  if(m_muonEfficiencySFTool){
    delete m_muonEfficiencySFTool;
    m_muonEfficiencySFTool = 0;
  }

  /// Muon Isolation Tool
  if(m_muonIsolationSFTool){
    delete m_muonIsolationSFTool;
    m_muonIsolationSFTool = 0;
  }

  /// Muon TTVA Efficiency Tool
  if(m_muonTTVAEfficiencySFTool){
    delete m_muonTTVAEfficiencySFTool;
    m_muonTTVAEfficiencySFTool = 0;
  }

  /// Muon Trigger Scale Factor Tool
  if(m_muonTriggerSFTool){
    delete m_muonTriggerSFTool;
    m_muonTriggerSFTool = 0;
  }

  /// Electron Efficiency Tool
  if(m_elecEfficiencySFTool_reco){
    delete m_elecEfficiencySFTool_reco;
    m_elecEfficiencySFTool_reco = 0;
  }

  if(m_elecEfficiencySFTool_id_Loose){
    delete m_elecEfficiencySFTool_id_Loose;
    m_elecEfficiencySFTool_id_Loose = 0;
  }

  if(m_elecEfficiencySFTool_id_Medium){
    delete m_elecEfficiencySFTool_id_Medium;
    m_elecEfficiencySFTool_id_Medium = 0;
  }

  if(m_elecEfficiencySFTool_id_Tight){
    delete m_elecEfficiencySFTool_id_Tight;
    m_elecEfficiencySFTool_id_Tight = 0;
  }

  if(m_elecEfficiencySFTool_iso_Loose){
    delete m_elecEfficiencySFTool_iso_Loose;
    m_elecEfficiencySFTool_iso_Loose = 0;
  }

  if(m_elecEfficiencySFTool_iso_Medium){
    delete m_elecEfficiencySFTool_iso_Medium;
    m_elecEfficiencySFTool_iso_Medium = 0;
  }

  if(m_elecEfficiencySFTool_iso_Tight){
    delete m_elecEfficiencySFTool_iso_Tight;
    m_elecEfficiencySFTool_iso_Tight = 0;
  }

  if(m_elecEfficiencySFTool_trigEff){
    delete m_elecEfficiencySFTool_trigEff;
    m_elecEfficiencySFTool_trigEff = 0;
  }

  if(m_elecEfficiencySFTool_trigSF){
    delete m_elecEfficiencySFTool_trigSF;
    m_elecEfficiencySFTool_trigSF = 0;
  }

  /// Jet JVT Efficiency Tool
  if(m_jvtefficiencyTool){
    delete m_jvtefficiencyTool;
    m_jvtefficiencyTool = 0;
  }

  /// Tau Efficiency Tool
  if(m_tauEffTool){
    delete m_tauEffTool;
    m_tauEffTool = 0;
  }

  /// MET Tools
  if(m_metSystTool){
    delete m_metSystTool;
    m_metSystTool = 0;
  }

  /// Isolation Correction Tool
  if(m_isoCorrTool){
    delete m_isoCorrTool;
    m_isoCorrTool = 0;
  }

  /// PileupReweighting Tool
  if(m_prwTool){
    delete m_prwTool;
    m_prwTool = 0;
  }

  /// Cutflow
  if(m_useBitsetCutflow && m_BitsetCutflow){
    delete m_BitsetCutflow;
    m_BitsetCutflow = 0;
  }


  // print out the number of Overlap removal
  Info("finalize()", "======================================================");
  Info("finalize()", "Number overlap electrons:    %i / %i", nOverlapElectrons, nInputElectrons);
  Info("finalize()", "Number overlap muons:    %i / %i", nOverlapMuons, nInputMuons);
  Info("finalize()", "Number overlap jets:    %i / %i", nOverlapJets, nInputJets);
  Info("finalize()", "Number overlap taus:    %i / %i", nOverlapTaus, nInputTaus);
  Info("finalize()", "Number overlap photons:    %i / %i", nOverlapPhotons, nInputPhotons);
  Info("finalize()", "======================================================");

  // print out the final number of clean events
  Info("finalize()", "Number of clean events = %i", m_numCleanEvents);

  // print out Cutflow
  Info("finalize()", "================================================");
  Info("finalize()", "================  VBF Cutflow  =================");
  for(int i=0; i<16 ; ++i) {
    int j = i+1;
    Info("finalize()", "Event cutflow (%i) = %i", j, m_eventCutflow[i]);
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MetTrigxAODAnalysis :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}



  EL::StatusCode MetTrigxAODAnalysis :: passMuonSelection(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(mu) = false;
    selectDec(mu) = false; // To select objects for Overlap removal
    dec_signal(mu) = false;

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double muPt = (mu.pt()) * 0.001; /// GeV
    //if ( muPt < 4. ) return EL::StatusCode::SUCCESS;

    // Muon Calibration
    if (!m_isData){
      if(m_muonCalibrationAndSmearingTool->applyCorrection(mu) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
        Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // MuonSelectionTool(Medium)
    if(!m_muonSelection->accept(mu)) return EL::StatusCode::SUCCESS;
    // MuonSelectionTool (Loose)
    //if(!m_loosemuonSelection->accept(mu)) return EL::StatusCode::SUCCESS;

    // Muon tranverse momentum cut
    if (muPt <= 10. ) return EL::StatusCode::SUCCESS; /// veto muon

    // Muon eta cut
    double muEta = mu.eta();
    if (fabs(muEta) >= 2.47) return EL::StatusCode::SUCCESS;

    //  if (mu.muonType()=xAOD::Muon_v1::Combined) return EL::StatusCode::SUCCESS;
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return EL::StatusCode::SUCCESS;


    // Baseline Muon
    dec_baseline(mu) = true;
    selectDec(mu) = true; // To select objects for Overlap removal


    // Muon pt cut
    if (muPt <= 25. ) return EL::StatusCode::SUCCESS;
    // Muon eta cut
    if (fabs(muEta) >= 2.4) return EL::StatusCode::SUCCESS;

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
      tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 3.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(mu)) return EL::StatusCode::SUCCESS;

    // Signal Muon
    dec_signal(mu) = true;

    return EL::StatusCode::SUCCESS;

  }



  EL::StatusCode MetTrigxAODAnalysis :: passMuonVBF(xAOD::Muon& mu,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(mu) = false;
    selectDec(mu) = false; // To select objects for Overlap removal

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passMuonSignal. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double muPt = (mu.pt()) * 0.001; /// GeV
    //if ( muPt < 4. ) return EL::StatusCode::SUCCESS;

    // Muon Calibration
    if (!m_isData){
      if(m_muonCalibrationAndSmearingTool->applyCorrection(mu) == CP::CorrectionCode::OutOfValidityRange){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
        Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
      }
    }

    // MuonSelectionTool (Loose)
    if(!m_loosemuonSelection->accept(mu)) return EL::StatusCode::SUCCESS;

    // Muon tranverse momentum
    if (muPt < m_muonPtCut ) return EL::StatusCode::SUCCESS;

    // Muon eta cut
    double muEta = mu.eta();
    if (fabs(muEta) > m_muonEtaCut) return EL::StatusCode::SUCCESS;

    // Combined (CB) or Segment-tagged (ST) muons (excluding Stand-alone (SA), Calorimeter-tagged (CaloTag) muons etc..)
    //if (!(mu.muonType() == xAOD::Muon::Combined || mu.muonType() == xAOD::Muon::SegmentTagged)) return EL::StatusCode::SUCCESS;
    if (mu.muonType() != xAOD::Muon_v1::Combined && mu.muonType() != xAOD::Muon_v1::SegmentTagged) return EL::StatusCode::SUCCESS;

    // d0 / z0 cuts applied
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle* tp;
    if (mu.muonType() == xAOD::Muon::SiliconAssociatedForwardMuon)
      tp = mu.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
    else
      tp = mu.primaryTrackParticle();
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 3.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( mu.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (muPt > 10. && muPt < 500. && !m_IsoToolVBF->accept(mu)) return EL::StatusCode::SUCCESS;

    dec_baseline(mu) = true;
    selectDec(mu) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }





  EL::StatusCode MetTrigxAODAnalysis :: passElectronSelection(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(elec) = false;
    selectDec(elec) = false; // To select objects for Overlap removal
    dec_signal(elec) = false;

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_identification

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double elecPt = (elec.pt()) * 0.001; /// GeV
    //if ( elecPt < 4. ) return EL::StatusCode::SUCCESS;

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return EL::StatusCode::SUCCESS;

    // "Please apply the identification to uncalibrated electron object. ID scale factors are to be applied to calibrated objects."
    // LH Electron (Medium)
    bool LHmediumSel = false;
    LHmediumSel = m_LHToolMedium2015->accept(elec);
    if (!LHmediumSel) return EL::StatusCode::SUCCESS;
    // LH Electron (Loose)
    //bool LHlooseSel = false;
    //LHlooseSel = m_LHToolLoose2015->accept(elec);
    //if (!LHlooseSel) return EL::StatusCode::SUCCESS;

    //Info("execute()", "  Selected electron pt from new Electron Container = %.2f GeV", (elec.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(elec) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original electron values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = elec.caloCluster()->eta();
    double Eta = elec.caloCluster()->etaBE(2);
    if ( fabs(Eta) >= 2.47 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

    // pT cut
    double elecPtCut = 10.0; /// GeV
    if (elecPt <= elecPtCut) return EL::StatusCode::SUCCESS; /// veto electron

    // Baseline Electron
    dec_baseline(elec) = true;
    selectDec(elec) = true; // To select objects for Overlap removal


    // pT cut
    if (elecPt <= 25.) return EL::StatusCode::SUCCESS; /// veto electron

    // d0 / z0 cuts applied
    // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle *tp = elec.trackParticle() ; //your input track particle from the electron
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 5.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( elec.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsolationSelectionTool->accept(elec)) return EL::StatusCode::SUCCESS;

    // Signal Electron
    dec_signal(elec) = true;

    return EL::StatusCode::SUCCESS;

  }



  EL::StatusCode MetTrigxAODAnalysis :: passElectronVBF(xAOD::Electron& elec,
      const xAOD::EventInfo* eventInfo, xAOD::Vertex* primVertex){

    dec_baseline(elec) = false;
    selectDec(elec) = false; // To select objects for Overlap removal

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passElectronSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // don't bother calibrating or computing WP
    double elecPt = (elec.pt()) * 0.001; /// GeV
    //if ( elecPt < 4. ) return EL::StatusCode::SUCCESS;

    // goodOQ(object quality cut) : Bad Electron Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !elec.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) return EL::StatusCode::SUCCESS;

    // "Please apply the identification to uncalibrated electron object. ID scale factors are to be applied to calibrated objects."
    // LH Electron identification
    //
    // LH Electron (Loose)
    bool LHlooseSel = false;
    LHlooseSel = m_LHToolLoose2015->accept(elec);
    if (!LHlooseSel) return EL::StatusCode::SUCCESS;
    /*
    // LH Electron (Medium)
    bool LHmediumSel = false;
    LHmediumSel = m_LHToolMedium2015->accept(elec);
    if (!LHmediumSel) return EL::StatusCode::SUCCESS;
    // LH Electron (Tight)
    bool LHtightSel = false;
    LHtightSel = m_LHToolTight2015->accept(elec);
    if (!LHtightSel) return EL::StatusCode::SUCCESS;
    */

    //Info("execute()", "  Selected electron pt from new Electron Container = %.2f GeV", (elec.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(elec) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original electron values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = elec.caloCluster()->eta();
    double Eta = elec.caloCluster()->etaBE(2);
    //if ( fabs(Eta) >= m_elecEtaCut || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;
    if ( fabs(Eta) > m_elecEtaCut ) return EL::StatusCode::SUCCESS;

    /// pT cut
    if (elecPt < m_elecPtCut ) return EL::StatusCode::SUCCESS; /// veto electron

    // d0 / z0 cuts applied
    // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio
    // d0 significance (Transverse impact parameter)
    const xAOD::TrackParticle *tp = elec.trackParticle() ; //your input track particle from the electron
    double d0sig = xAOD::TrackingHelpers::d0significance( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if (fabs(d0sig) > 5.0) return EL::StatusCode::SUCCESS;
    // zo cut
    float z0sintheta = 1e8;
    //if (primVertex) z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( elec.p4().Theta() );
    z0sintheta = ( tp->z0() + tp->vz() - primVertex->z() ) * TMath::Sin( tp->theta() );
    if (fabs(z0sintheta) > 0.5) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    if (!m_IsoToolVBF->accept(elec)) return EL::StatusCode::SUCCESS;

    dec_baseline(elec) = true;
    selectDec(elec) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }





  EL::StatusCode MetTrigxAODAnalysis :: passPhotonSelection(xAOD::Photon& phot,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(phot) = false;
    selectDec(phot) = false; // To select objects for Overlap removal

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passPhotonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Photon author cuts
    if ( !(phot.author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) )
      return EL::StatusCode::SUCCESS;


    //Info("execute()", "  Selected photon pt from new Photon Container = %.2f GeV", (phot.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = phot.caloCluster()->eta();
    double Eta = phot.caloCluster()->etaBE(2);
    if ( fabs(Eta) >= 2.37 || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;

    // pT cut
    double photPt = (phot.pt()) * 0.001; /// GeV
    double photPtCut = 10.0; /// GeV
    if (photPt <= photPtCut) return EL::StatusCode::SUCCESS; /// veto photon

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return EL::StatusCode::SUCCESS;

    // MC fudge tool
    if (!m_isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    //if (!m_IsolationSelectionTool->accept(phot)) return EL::StatusCode::SUCCESS;


    dec_baseline(phot) = true;
    selectDec(phot) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }


  EL::StatusCode MetTrigxAODAnalysis :: passPhotonVBF(xAOD::Photon& phot,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(phot) = false;
    selectDec(phot) = false; // To select objects for Overlap removal

    // According to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passPhotonSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Photon author cuts
    if ( !(phot.author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) )
      return EL::StatusCode::SUCCESS;


    //Info("execute()", "  Selected photon pt from new Photon Container = %.2f GeV", (phot.pt() * 0.001));

    // Calibration
    if(m_egammaCalibrationAndSmearingTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
      // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
      // If OutOfValidityRange is returned no modification is made and the original photon values are taken.
      Error("execute()", "EgammaCalibrationAndSmearingTool returns Error CorrectionCode");
    }

    // Eta cut
    //double Eta = phot.caloCluster()->eta();
    double Eta = phot.caloCluster()->etaBE(2);
    //if ( fabs(Eta) >= m_photEtaCut || (fabs(Eta) >= 1.37 && fabs(Eta) <= 1.52)) return EL::StatusCode::SUCCESS;
    if ( fabs(Eta) > m_photEtaCut ) return EL::StatusCode::SUCCESS;

    // pT cut
    double photPt = (phot.pt()) * 0.001; /// GeV
    if (photPt < m_photPtCut) return EL::StatusCode::SUCCESS; /// veto photon

    // goodOQ(object quality cut) : Bad photon Cluster
    // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EGammaIdentificationRun2#Object_quality_cut
    if( !phot.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON) ) return EL::StatusCode::SUCCESS;

    // MC fudge tool
    if (!m_isData){
      if(m_electronPhotonShowerShapeFudgeTool->applyCorrection(phot) == CP::CorrectionCode::Error){ // apply correction and check return code
        Error("execute()", "ElectronPhotonShowerShapeFudgeTool returns Error CorrectionCode");
      }
    }

    // Recomputing the photon ID flags
    if (!m_photonTightIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;
    //if (!m_photonLooseIsEMSelector->accept(phot)) return EL::StatusCode::SUCCESS;

    // Isolation requirement
    //if (!m_IsoToolVBF->accept(phot)) return EL::StatusCode::SUCCESS;


    dec_baseline(phot) = true;
    selectDec(phot) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }



  EL::StatusCode MetTrigxAODAnalysis :: passTauSelection(xAOD::TauJet& tau,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(tau) = false;
    selectDec(tau) = false; // To select objects for Overlap removal

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // Tau Smearing (for MC)
    if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !m_isData){ // it's MC!
      if(m_tauSmearingTool->applyCorrection(tau) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original tau values are taken.
        Error("execute()", "TauSmearingTool returns Error CorrectionCode");
      }
    }

    //Info("execute()", "  original tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    // TauSelectionTool (Medium)
    if(!m_tauSelTool->accept(tau)) return EL::StatusCode::SUCCESS;
    // TauSelectionTool (Loose for VBF)
    //if(!m_tauSelToolVBF->accept(tau)) return EL::StatusCode::SUCCESS;

    //Info("execute()", "  Selected tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    dec_baseline(tau) = true;
    selectDec(tau) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }


  EL::StatusCode MetTrigxAODAnalysis :: passTauVBF(xAOD::TauJet& tau,
      const xAOD::EventInfo* eventInfo){

    dec_baseline(tau) = false;
    selectDec(tau) = false; // To select objects for Overlap removal

    // According to https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/README.rst

    // Event information
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
      Error("execute()", "Failed to retrieve event info collection in passTauSelection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // TauOverlappingElectronLLHDecorator
    m_tauOverlappingElectronLLHDecorator->decorate(tau);

    // Tau Smearing (for MC)
    if( fabs(tau.eta()) <= 2.5 && tau.nTracks() > 0 && !m_isData){ // it's MC!
      if(m_tauSmearingTool->applyCorrection(tau) == CP::CorrectionCode::Error){ // apply correction and check return code
        // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
        // If OutOfValidityRange is returned no modification is made and the original tau values are taken.
        Error("execute()", "TauSmearingTool returns Error CorrectionCode");
      }
    }

    //Info("execute()", "  original tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    // TauSelectionTool (Loose for VBF)
    if(!m_tauSelToolVBF->accept(tau)) return EL::StatusCode::SUCCESS;

    //Info("execute()", "  Selected tau pt from new Tau Container = %.2f GeV", (tau.pt() * 0.001));

    dec_baseline(tau) = true;
    selectDec(tau) = true; // To select objects for Overlap removal

    return EL::StatusCode::SUCCESS;

  }




  bool MetTrigxAODAnalysis :: IsBadJet(xAOD::Jet& jet) {

    if (m_doORtool && overlapAcc(jet)) return false;

    double jetPt = (jet.pt()) * 0.001; /// GeV

    //Info("execute()", "  corrected jet pt in IsBadJet function = %.2f GeV", jetPt );
    //Info("execute()", "  updated jet jvt in IsBadJet function = %.2f ", cacc_jvt(jet) );

    // Pile-up
    if ( cacc_jvt(jet) < 0.59 && fabs(jet.eta()) < 2.4 && jetPt < 50.0 ) return false;

    // pT cut
    if ( jetPt < m_jetPtCut || fabs(jet.eta()) > m_jetEtaCut) return false; 

    // Jet Cleaning Tool
    dec_bad(jet) = !m_jetCleaningLoose->accept( jet );

    return dec_bad(jet);

  }


  bool MetTrigxAODAnalysis :: IsSignalJet(xAOD::Jet& jet) {

    if (m_doORtool && overlapAcc(jet)) return false;
    if (!dec_baseline(jet)) return false;

    double jetPt = (jet.pt()) * 0.001; /// GeV

    // pT, eta cut
    if ( jetPt < m_jetPtCut || fabs(jet.eta()) > m_jetEtaCut) return false;

    bool isgoodjet = !dec_bad(jet) && (cacc_jvt(jet) > 0.59 || fabs(jet.eta()) > 2.4 || jetPt > 50.0);

    dec_signal(jet) = isgoodjet;

    return isgoodjet;

  }


  float MetTrigxAODAnalysis :: GetGoodMuonSF(xAOD::Muon& mu,
      const bool recoSF, const bool isoSF, const bool ttvaSF) {

    float sf(1.);

    if (recoSF) {
      float sf_reco(1.);
      if (m_muonEfficiencySFTool->getEfficiencyScaleFactor( mu, sf_reco ) == CP::CorrectionCode::Error) {
        Error("execute()", " GetGoodMuonSF: Reco getEfficiencyScaleFactor returns Error CorrectionCode");
      }
      //Info("execute()", "  GetGoodMuonSF: sf_reco = %.5f ", sf_reco );
      sf *= sf_reco;
    }

    if (isoSF) {
      float sf_iso(1.);
      if (m_muonIsolationSFTool->getEfficiencyScaleFactor( mu, sf_iso ) == CP::CorrectionCode::Error) {
        Error("execute()", " GetGoodMuonSF: Iso getEfficiencyScaleFactor returns Error CorrectionCode");
      }
      //Info("execute()", "  GetGoodMuonSF: sf_iso = %.5f ", sf_iso );
      sf *= sf_iso;
    }

    if (ttvaSF) {
      float sf_TTVA(1.);
      if (m_muonTTVAEfficiencySFTool->getEfficiencyScaleFactor( mu, sf_TTVA ) == CP::CorrectionCode::Error) {
        Error("execute()", " GetGoodMuonSF: TTVA getEfficiencyScaleFactor returns Error CorrectionCode");
      }
      //Info("execute()", "  GetGoodMuonSF: sf_TTVA = %.5f ", sf_TTVA );
      sf *= sf_TTVA;
    }

    //Info("execute()", "  GetGoodMuonSF: Good Muon SF = %.5f ", sf );
    dec_scalefactor(mu) = sf;
    return sf;

  }

  double MetTrigxAODAnalysis :: GetTotalMuonSF(xAOD::MuonContainer& muons,
      bool recoSF, bool isoSF, bool ttvaSF) {

    double sf(1.);

    for (const auto& muon : muons) {
      double muPt = muon->pt() * 0.001; /// GeV
      if (muPt < m_isoMuonPtMin || muPt > m_isoMuonPtMax) {
        isoSF = false;
      }
      //Info("execute()", "  GetTotalMuonSF: Muon pT = %.2f GeV ", muPt );
      sf *= GetGoodMuonSF(*muon, recoSF, isoSF, ttvaSF);
    }

    //Info("execute()", "  GetTotalMuonSF: Total Muon SF = %.5f ", sf );
    return sf;

  }

  int MetTrigxAODAnalysis :: NumIsoTracks(const xAOD::TrackParticleContainer* inTracks,
      xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
    //
    //  Fill track objects with information about isolated tracks. For being isolated, there should be no track
    //  above 3 GeV satisfying quality requirement that are within a cone of 0.4 around the probed track.
    //
    //============================================================================================================

    // Integer to return with this function
    // ------------------------------------

    int NisoTrack = 0;


    // Loop over tracks in the event
    // -----------------------------

    for( auto trk_itr : *inTracks ){

      int NCloseby = 0;

      // Apply quality cuts on the track with Pt>10 GeV
      // ----------------------------------------------
      float pt   = (trk_itr->pt()) * 0.001; /// GeV
      float eta  = trk_itr->eta();
      //float phi  = trk_itr->phi();
      float d0   = trk_itr->d0();
      float z0   = (trk_itr->z0() + trk_itr->vz() - primVertex->z());
      int Ndof    = trk_itr->numberDoF();
      if (Ndof==0) Ndof = 1;

      float Chi2 = trk_itr->chiSquared()/Ndof;
      uint8_t nSCT      = -1;
      uint8_t nPix      = -1;
      if(!trk_itr->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
      if(!trk_itr->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
      uint8_t NHits   = nSCT + nPix;

      if(pt < Pt_High || fabs(eta) > 2.5 || fabs(z0) >= 2.0 || fabs(d0) >= 1.0 || Chi2 >= 3.0 || NHits < 5) continue;

      // Loop over *other* tracks and apply quality cuts too
      // ---------------------------------------------------

      for( auto trk2_itr : *inTracks ){

        if ( trk2_itr == trk_itr ) continue;

        float pt_2   = (trk2_itr->pt()) * 0.001; /// GeV
        float eta_2  = trk2_itr->eta();
        //float phi_2  = trk2_itr->phi();
        float d0_2   = trk2_itr->d0();
        float z0_2   = (trk2_itr->z0() + trk2_itr->vz() - primVertex->z());
        int Ndof_2    = trk2_itr->numberDoF();
        if (Ndof_2==0) Ndof_2 = 1;

        float Chi2_2 = trk2_itr->chiSquared()/Ndof_2;
        uint8_t nSCT_2      = -1;
        uint8_t nPix_2      = -1;
        if(!trk2_itr->summaryValue(nPix_2,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
        if(!trk2_itr->summaryValue(nSCT_2,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
        uint8_t NHits_2   = nSCT_2 + nPix_2;

        if(pt_2 < Pt_Low || fabs(eta_2) > 2.5 || fabs(z0_2) >= 2.0 || fabs(d0_2) >= 1.0 || Chi2_2 >= 3.0 || NHits_2 < 5) continue;

        float dR = trk2_itr->p4().DeltaR(trk_itr->p4());
        //Info("execute()", "  dR(track1 with pt=%f, track2 with pt=%f) = %f ", pt, pt_2, dR);

        // Count the number of track in a 0.4 cone around the probed track
        // ---------------------------------------------------------------

        if(dR < 0.4) NCloseby++;
      } // end 2nd loop

      // Fill the vertex collection
      // --------------------------

      // Note: The object properties are: E, pT, eta, phi, d0, d0_err, d0_vx, d0_err_vx, charge, n_SCTHits, n_PixelHits

      if (NCloseby < 1) NisoTrack++;

    } // end 1st loop

    return NisoTrack;
  }


  int MetTrigxAODAnalysis :: NumMuonIsoTrack(xAOD::MuonContainer* muons, const xAOD::TrackParticleContainer* inTracks,
      xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
    //
    //  Apply the same criteria as in the SetIsoTracks function to determine how many muons are isolated
    //  according to this definition.
    //
    //============================================================================================================

    // Integer to return with this function
    // ------------------------------------

    int NisoMuon = 0;


    // Loop over muon objects obtained in the event
    // --------------------------------------------

    for (const auto& muon_itr : *muons) { // C++11 shortcut

      int NCloseby = 0;

      // Apply muon cuts on the baseline muon with pT>10GeV
      float muon_pt = (muon_itr->pt()) * 0.001; /// GeV

      //if (!dec_baseline(muon_itr) || !m_IsolationSelectionTool->accept(muon_itr)) continue; 
      if (!dec_baseline(*muon_itr) || muon_pt < Pt_High) continue; 

      // Loop over *other* tracks and apply quality cuts too
      // ---------------------------------------------------

      for( auto trk_itr : *inTracks ){

        // Apply quality cuts on the track with Pt>3 GeV
        // ----------------------------------------------
        float pt   = (trk_itr->pt()) * 0.001; /// GeV
        float eta  = trk_itr->eta();
        //float phi  = trk_itr->phi();
        float d0   = trk_itr->d0();
        float z0   = (trk_itr->z0() + trk_itr->vz() - primVertex->z());
        int Ndof    = trk_itr->numberDoF();
        if (Ndof==0) Ndof = 1;

        float Chi2 = trk_itr->chiSquared()/Ndof;
        uint8_t nSCT      = -1;
        uint8_t nPix      = -1;
        if(!trk_itr->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
        if(!trk_itr->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
        uint8_t NHits   = nSCT + nPix;

        if(pt < Pt_Low || fabs(eta) > 2.5 || fabs(z0) >= 2.0 || fabs(d0) >= 1.0 || Chi2 >= 3.0 || NHits < 5) continue;


        float dR = trk_itr->p4().DeltaR(muon_itr->p4());
        //Info("execute()", "  dR(Muon with pt=%f, track with pt=%f) = %f ", muon_pt, pt, dR);

        // Count the number of track in a 0.4 cone around the probed track
        // ---------------------------------------------------------------

        if(dR < 0.4) NCloseby++;
      } // end 2nd loop

      // Count the number of isolated muon in the event following track iso criteria
      // ---------------------------------------------------------------------------
      //
      // Note 1: We assume here that the track in the collection which is the closest to the muon is actually the
      //         one that come from the muon. This will happen for each muon, so the muon will be isolated if there
      //         other tracks except this clostest one which is the muon itself (hence NCloseby<=1). It is also
      //         possible that no tracks are found in the cone because of the Chi2 cut on the muon might not exactly
      //         match the one required for the closeby tracks.
      //
      // Note 2: Typically, there will be one Iso muon per W events and two for Z.

      if (NCloseby <= 1) NisoMuon++;

    } // end 1st loop

    return NisoMuon;
  }


  int MetTrigxAODAnalysis :: NumElecIsoTrack(xAOD::ElectronContainer* electrons, const xAOD::TrackParticleContainer* inTracks,
      xAOD::Vertex* primVertex, float Pt_Low, float Pt_High) {
    //
    //  Apply the same criteria as in the SetIsoTracks function to determine how many electrons are isolated
    //  according to this definition.
    //
    //============================================================================================================

    // Integer to return with this function
    // ------------------------------------

    int NisoElec = 0;


    // Loop over electron objects obtained in the event
    // --------------------------------------------

    for (const auto& elec_itr : *electrons) { // C++11 shortcut

      int NCloseby = 0;

      // Apply electron cuts on the baseline electron with pT>10GeV
      float elec_pt = (elec_itr->pt()) * 0.001; /// GeV

      //if (!dec_baseline(elec_itr) || !m_IsolationSelectionTool->accept(elec_itr)) continue; 
      if (!dec_baseline(*elec_itr) || elec_pt < Pt_High) continue; 

      // Loop over *other* tracks and apply quality cuts too
      // ---------------------------------------------------

      for( auto trk_itr : *inTracks ){

        // Apply quality cuts on the track with Pt>3 GeV
        // ----------------------------------------------
        float pt   = (trk_itr->pt()) * 0.001; /// GeV
        float eta  = trk_itr->eta();
        //float phi  = trk_itr->phi();
        float d0   = trk_itr->d0();
        float z0   = (trk_itr->z0() + trk_itr->vz() - primVertex->z());
        int Ndof    = trk_itr->numberDoF();
        if (Ndof==0) Ndof = 1;

        float Chi2 = trk_itr->chiSquared()/Ndof;
        uint8_t nSCT      = -1;
        uint8_t nPix      = -1;
        if(!trk_itr->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("PassCuts()", "Pix hits not filled");
        if(!trk_itr->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("PassCuts()", "SCT hits not filled");
        uint8_t NHits   = nSCT + nPix;

        if(pt < Pt_Low || fabs(eta) > 2.5 || fabs(z0) >= 2.0 || fabs(d0) >= 1.0 || Chi2 >= 3.0 || NHits < 5) continue;


        float dR = trk_itr->p4().DeltaR(elec_itr->p4());
        //Info("execute()", "  dR(Electron with pt=%f, track with pt=%f) = %f ", elec_pt, pt, dR);

        // Count the number of track in a 0.4 cone around the probed track
        // ---------------------------------------------------------------

        if(dR < 0.4) NCloseby++;
      } // end 2nd loop

      // Count the number of isolated elec in the event following track iso criteria
      // ---------------------------------------------------------------------------

      if (NCloseby <= 1) NisoElec++;

    } // end 1st loop

    return NisoElec;
  }



  float MetTrigxAODAnalysis :: DeltaPhi(float phi1, float phi2) {

    float dPhi = TMath::Abs(phi1 - phi2);

    if(dPhi > TMath::Pi())
      dPhi = TMath::TwoPi() - dPhi;

    return dPhi;

  }



  float MetTrigxAODAnalysis :: DeltaR(float eta1, float eta2, float phi1, float phi2) {

    float dEta = eta1 - eta2;
    float dPhi = DeltaPhi(phi1,phi2);

    return TMath::Sqrt(dEta*dEta + dPhi*dPhi);

  }
