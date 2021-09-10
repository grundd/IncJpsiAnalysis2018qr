
// c++ headers
#include <iostream>
#include <fstream>

// root headers
#include <TMath.h>


double nextPt(double f, double b, double ti, double tf, double t1)
{

  // exponatial
  double ei = TMath::Exp(-b*ti);
  double ef = TMath::Exp(-b*tf);
  double e1 = TMath::Exp(-b*t1);

  // integral in the full range
  double iT = ei-ef;

  // new t
  double tn = -TMath::Log(e1-f*iT)/b;
  
  //  this part is to test
  bool testing = kFALSE;
  if (testing) {
    // integral over new bin
    double iB = e1-TMath::Exp(-b*tn);
    
    std::cout << " pTn = " << TMath::Sqrt(tn) << std::endl;
    std::cout << " iT = " << iT << std::endl;
    std::cout << " iB = " << iB << std::endl;
    std::cout << " iT/iB = " << (iT/iB) << std::endl;
  }  // end of testing

  return TMath::Sqrt(tn);
}

void ptBinning(int nBins, float b = 4,
	       float pTi = 0.2, float pTf = 1.0)
{
  // Mandeltan t
  double ti =  pTi*pTi; // initial t
  double tf =  pTf*pTf; // final t
  
  // initialise
  double t1 = ti;

  // fraction ocupied by each bin
  double f = 1.0/nBins;

  // find each new bin boarder
  for(int i=nBins;i>0;i--) {
    double ptn = nextPt(f,b,ti,tf,t1);
    std::cout << " bin " << (nBins-i+1)
	      <<" = ("<<TMath::Sqrt(t1)
	      <<","<<ptn<<")" <<std::endl;
    // prepare for next bin
    t1 = ptn*ptn;
  }
}
