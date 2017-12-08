#ifndef GEOM_H
#define GEOM_H


#include "TVector3.h"
#include "TMath.h"

Double_t dist3D(TVector3 p, TVector3 u, TVector3 q, TVector3 v, 
		TVector3 &pc1, TVector3 &pc2){
  //line equation:
  // l1 = p + u * s
  // l2 = q + v * t
  //s, t are the "time" parameters

  TVector3 w = p-q;
  double   a = u*u;         
  double   b = u*v;
  double   c = v*v;         
  double   d = u*w;
  double   e = v*w;
  double   D = a*c - b*b;
  double   sc, tc;
  
  // compute the line parameters of the two closest points
  if (D < 1E-6) {          // if the lines are almost parallel
    sc = 0.0;
    tc = (b>c ? d/b : e/c);    // then use the largest denominator
  }
  else {
    sc = (b*e - c*d) / D;
    tc = (a*e - b*d) / D;
  }
  
  // get the difference of the two closest points
  // dist = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)

  pc1 = p + (sc * u);
  pc2 = q + (tc * v);
  
  return (pc1-pc2).Mag();
}

void getTrack2D(Double_t b,Double_t phi, Double_t rCell,
		TVector3&p1, TVector3&p2){

  Double_t x0=b*sin(phi); //cos(phi-TMath::PiOver2());
  Double_t y0=b*(-cos(phi)); //=sin(phi-TMath::PiOver2());

  Double_t x1=cos(phi); 
  Double_t y1=sin(phi); 

  // track equation:
  // x=x0+t*x1
  // y=y0+t*y1 -->
  // y=y0 + (x-x0)*(y1/x1)

  Double_t m = y1/x1;
  Double_t n = y0 - x0*m;

    // track equation
    // y=m*x+n
    // or x=y/m-n/m
    Double_t y_a = rCell*m + n;
    Double_t y_b = -rCell*m + n;
    Double_t x_a = rCell/m - n/m;
    Double_t x_b = -rCell/m - n/m;

    p1.SetXYZ(0,0,0);
    p2.SetXYZ(0,0,0);
    TVector3 *p=&p1;

    if(y_a>=-rCell && y_a<=rCell){ 
      p->SetXYZ(rCell,y_a,0); 
      p=&p2;
    }
    if(y_b>=-rCell && y_b<=rCell){ 
      p->SetXYZ(-rCell,y_b,0);
      p=&p2; 
    }
    if(x_a>-rCell && x_a<rCell){ 
      p->SetXYZ(x_a,rCell,0); 
      p=&p2;
    }
    if(x_b>-rCell && x_b<rCell){ 
      p->SetXYZ(x_b,-rCell,0); 
      p=&p2;
    }

}

#endif
