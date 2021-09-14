void RunLumiAnalysis()
{
  gROOT->ProcessLine(".include $ALICE_ROOT/include");

  gROOT->ProcessLine(".x LuminosityCalculationStandAlone.C+g");

  Printf("Lumi calculation finished.");
}
