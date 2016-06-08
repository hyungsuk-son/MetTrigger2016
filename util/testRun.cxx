#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>
#include "SampleHandler/ScanDir.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

#include <MetTriggerPackage/MetTrigxAODAnalysis.h>

int main( int argc, char* argv[] ) {

  // Take the submit directory from the input if provided:
  std::string submitDir = "submitDir";
  if( argc > 1 ) submitDir = argv[ 1 ];

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // Construct the samples to run on:
  SH::SampleHandler sh;

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  //const char* inputFilePath = gSystem->ExpandPathName ("$ALRB_TutorialData/r6630/");
  //SH::ScanDir().sampleDepth(1).samplePattern("AOD.05352803._000031.pool.root.1").scan(sh, inputFilePath);
  const char* inputFilePath = gSystem->ExpandPathName ("/cluster/tufts/atlas08/DATA/Run2");
  //SH::ScanDir().filePattern("AOD.06828599._000351.pool.root.1").scan(sh,inputFilePath); // Data15
  SH::ScanDir().filePattern("AOD.08517453._*.pool.root.1").scan(sh,inputFilePath); // Data16
  //SH::ScanDir().filePattern("DAOD_EXOT5.07992543._000007.pool.root.1").scan(sh,inputFilePath); // MC Zvv
  //SH::ScanDir().filePattern("DAOD_EXOT5.07992459._000012.pool.root.1").scan(sh,inputFilePath); // MC Zmumu


  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  //job.options()->setDouble (EL::Job::optMaxEvents, 500); // for testing
/*
  // For ntuple
  // define an output and an ntuple associated to that output (For ntuple)
  EL::OutputStream output  ("myOutput");
  job.outputAdd (output);
  EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
  job.algsAdd (ntuple);
  //---------------------------------------------------------------------
*/
  // Add our analysis to the job:
  MetTrigxAODAnalysis* alg = new MetTrigxAODAnalysis();
  job.algsAdd( alg );
/*
  // For ntuple
  // Let your algorithm know the name of the output file (For ntuple)
  alg->outputName = "myOutput"; // give the name of the output to our algorithm
  //----------------------------------------------------------------------------
*/
  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  driver.submit( job, submitDir );

  return 0;
}
