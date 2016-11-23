#ifndef MYMATH
#define MYMATH

namespace MyMath {

  double DeltaPhi(double phi1, double phi2) 
  {
    
    double dphi = phi2 - phi1; 

    if ( dphi > M_PI ) {
      dphi -= 2.0*M_PI; 
    } 
    else if ( dphi <= -M_PI ) {
      dphi += 2.0*M_PI;
    }

    return dphi;

  }

}

#endif
