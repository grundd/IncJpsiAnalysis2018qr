void RunLumiAnalysis()
{
    gROOT->ProcessLine(".L LuminosityCalculationStandAlone.C+");

    gROOT->ProcessLine(".x LuminosityCalculationStandAlone.C\(0\)");

    gROOT->ProcessLine(".x LuminosityCalculationStandAlone.C\(1\)");

    Printf("Lumi calculation finished.");
}
