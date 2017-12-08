#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TDatime.h"
#include "TVector2.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "DriftLineRKF.hh"
#include "TROOT.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRandom3.h"
#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "SolidTube.hh"
#include "ComponentAnalyticField.hh"
#include "GeometrySimple.hh"
#include "ViewCell.hh"
#include "ViewGeometry.hh"
#include "ViewSignal.hh"
#include "ViewField.hh"
#include "ViewDrift.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "ViewMedium.hh"
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include <iostream>
#include "geom.h"

using namespace Garfield;
using namespace std;

int main(int argc,char *argv[]){

  TString gasFile = "../data/gas/He_95_Ar_1_Eth_4.gas";
  //TString gasFile = "../data/gas/Ar_90_CO2_10.gas";
  int Ntrk=2;
  Bool_t view = false;
  // Cell radius [cm]
  float rCell = 2.2;
  float gapWire =  0.911;
  // Wire radius [cm]
  const double rSWire = 0.002;
  const double rFWire = 0.0025;
  // Half-length of the tube [cm]
  const double lTube = 360.;
  // Voltages
  const double vWire = 3000.;

  float energy = 2995.*1000./2.; // in eV
  TString particle = "e-";

  for(Int_t iar = 1; iar<argc; iar++){

    TString arg = argv[iar];
    
    cout << arg << endl;
    if(arg.BeginsWith("--gas=")){
      arg.ReplaceAll("--gas=","");
      gasFile = arg;
    } else if(arg.BeginsWith("--N=")){
      arg.ReplaceAll("--N=","");
      cout << arg << endl;
      Ntrk = arg.Atoi();
    } else if(arg.BeginsWith("--rCell=")){
       arg.ReplaceAll("--rCell=","");
       float rCell2 = arg.Atof();
       if(rCell2<=0) {
	 cout << " rCell = " << rCell << endl;
	 return 1;
       }
       float f = rCell2/rCell;
       gapWire = f*gapWire;
       rCell = rCell2; // = f*rCell;
    } else if(arg.BeginsWith("--E=")){
      arg.ReplaceAll("--E=","");
      energy = arg.Atof()*1E6;
      if(energy<=0) {
	cout << " energy = " << energy << endl;
	  return 1;
      }
    } else if(arg.BeginsWith("--part=")){
      arg.ReplaceAll("--part=","");
      particle = arg;
    } else if(arg.BeginsWith("--view")){
      view = true;
    } else if(arg.BeginsWith("-h") | arg.BeginsWith("--help")){
	cout << "clugen options:" << endl
	     << "\t--gas=<gas_file> " << endl
	     << "\t--part=<particle> " << endl
	     << "\t--E=<particle energy in MeV> " << endl
	     << "\t--N=<number_of_events>" << endl
	     << "\t--rCell=<cell_radius>" << endl
	     << "\t--view \t\t [enable plots]" << endl
	     << "\t--help|-h \t [this screen]" << endl
	     // << "\t" << endl
	     << "" << endl;
	return 1;
    } else {
      cout << arg << " argument not recognised (use -h for help)" << endl;
      return 1;
    }
  }

  TDatime date;
  gRandom->SetSeed(date.GetTime());
  TApplication *app;
  if(view) app = new TApplication("app",0,NULL);

  const double sqrt2 = sqrt(2);
  //setup gas:
  MediumMagboltz* gas = new MediumMagboltz();
  gas->LoadGasFile(gasFile.Data());
  
  /////SETUP THE GEOMETRY
  ComponentAnalyticField* cmp = new ComponentAnalyticField();
  GeometrySimple* geo = new GeometrySimple();

  TString title = Form("p=(%s)-E=%0.0fMeV-v=%0.0f-rCell=%0.2f-%s",particle.Data(),energy/1E6,vWire,rCell,gSystem->BaseName(gasFile));

  cout << title << endl;

  //Make a large tube of gas including all the chamber
  SolidTube *tube = new SolidTube(0.,0.,0.,rSWire,2*rCell,lTube);
  geo->AddSolid(tube,gas);

  cmp->SetGeometry(geo);
  //Add sense wires
  cmp->AddWire(0,0,rSWire,vWire,"s"); 
  //Add field wires
  cmp->AddWire(-rCell,0,rFWire,0.,"f");
  cmp->AddWire(-rCell,-gapWire,rFWire,0.,"f");
  cmp->AddWire(-rCell,gapWire,rFWire,0.,"f");

  cmp->AddWire(0,rCell,rFWire,0.,"f");
  cmp->AddWire(-gapWire,rCell,rFWire,0.,"f");
  cmp->AddWire(gapWire,rCell,rFWire,0.,"f");

  cmp->SetPeriodicityX(rCell*2);
  cmp->SetPeriodicityY(rCell*2);

  cmp->AddReadout("s");
  cmp->AddReadout("f");

  //Define a sensor
  Sensor* sensor = new Sensor();
  sensor->AddComponent(cmp);
  sensor->AddElectrode(cmp, "s");
  sensor->SetTimeWindow(0., 0.1, 1000);
  sensor->SetArea(-2*rCell,-2*rCell,-lTube,2*rCell,2*rCell,lTube);
  const double tMin = 0.;
  const double tMax = 25000.;
  const double tStep = 0.02;
  const int nTimeBins = int((tMax - tMin) / tStep);
  sensor->SetTimeWindow(0., tStep, nTimeBins);

  //Tracking and avalanches
  AvalancheMC* aval = new AvalancheMC(); 
  aval->SetSensor(sensor); 
  aval->EnableDiffusion();
  aval->EnableSignalCalculation();
	
	
  TrackHeed *Tracking = new TrackHeed();
  Tracking->SetParticle(particle.Data());
  Tracking->SetEnergy(energy);
  Tracking->SetSensor(sensor);

  ViewDrift* view_drift = 0;
  if(view){
    view_drift = new ViewDrift(); 
    view_drift->SetArea(-2*rCell,-2*rCell,-lTube,2*rCell,2*rCell,lTube);
    Tracking->EnablePlotting(view_drift);
  }

#define  MAXCLS 50

  Int_t nclu, ne[MAXCLS];
  Int_t nelec[MAXCLS];
  Double_t cls[MAXCLS][3], tcls[MAXCLS], ecls[MAXCLS], dhit[MAXCLS], thit[MAXCLS];
  Double_t orig[3],  dir[3], p1[3], p2[3];
  Double_t b, theta, phi;
  Double_t integral[MAXCLS];
  Double_t gain[MAXCLS];

  TFile *sfile = new TFile(title+Form("-%dclusters.root",Ntrk),"RECREATE");
  cout << sfile->GetName() << endl;
  TTree *tree = new TTree("clusters","clusters");
  tree->Branch("b",&b,"b/D");
  tree->Branch("theta",&theta,"theta/D");
  tree->Branch("phi",&phi,"phi/D");
  tree->Branch("nclu",&nclu,"nclu/I");
  tree->Branch("integral",&integral,"integral[nclu]/D");
  tree->Branch("gain",&gain,"gain[nclu]/D");
	
  tree->Branch("ne",ne,"ne[nclu]/I");
  tree->Branch("nelec",nelec,"nelec[nclu]/I");
  tree->Branch("cls",cls,"cls[nclu][3]/D");
  tree->Branch("tcls",tcls,"tcls[nclu]/D");
  tree->Branch("ecls",ecls,"ecls[nclu]/D");
  tree->Branch("dhit",dhit,"dhit[nclu]/D");
  tree->Branch("thit",thit,"thit[nclu]/D");

  tree->Branch("orig",orig,"orig[3]/D");
  tree->Branch("dir",dir,"dir[3]/D");
  tree->Branch("p1",p1,"p1[3]/D");
  tree->Branch("p2",p2,"p2[3]/D");


  for(Int_t itrk=0;itrk<Ntrk;itrk++){
    // cout << endl << "event " << itrk << endl;
    if(itrk%10==0) cout << itrk << endl;
    TVector3 vorig, vdir;
    TVector3 vp1,vp2,vzero(0,0,0);
	  ////
	  b = gRandom->Uniform(0,sqrt2*rCell);
	  phi=gRandom->Uniform(0,TMath::TwoPi());
	  theta=TMath::PiOver2();
	  
	  getTrack2D(b,phi,rCell,vp1,vp2);
	  if( (vp1==vzero) | (vp2==vzero)){ itrk--; continue;}
	  
	  vorig=vp1;
	  vdir=vp2-vp1;
	  vdir.SetTheta(theta);
	  
	  vorig.GetXYZ(orig);
	  vdir.GetXYZ(dir);
	  vp1.GetXYZ(p1);
	  vp2.GetXYZ(p2);
	  
	 
	  
	  // TVector3 vwire(0,0,1),va,vb;
	  // cout << "b = " << b;
	  // 	 << " =?= " << dist3D(vzero,vwire,vorig,vdir,va,vb)
	  // 	 << " phi = " << phi
	  // 	 << " theta = " << theta
	  // 	 << endl;
	  // cout << "vorig "; vorig.Print();
	  // cout << "vdir ";vdir.Print();
	  // cout << "vp1 ";vp1.Print();
	  // cout << "vp2 ";vp2.Print();
	  // cout << "va ";va.Print();
	  // cout << "vb ";vb.Print();
	  
	  ////
    Tracking->NewTrack(vorig.X(),vorig.Y(),vorig.Z(), //Position
		       0.,                   //Time
		       vdir.X(),vdir.Y(),vdir.Z());  //Direction

    Double_t extra_;
	  
    nclu = 0;
    while(Tracking->GetCluster(cls[nclu][0],cls[nclu][1],cls[nclu][2],
			       tcls[nclu],
			       ne[nclu],ecls[nclu],extra_)){
       cout << nclu << " - " << cls[nclu][0] << ", " << cls[nclu][1]<< ", " << cls[nclu][2]<< ", " << tcls[nclu]<< ", " << ne[nclu]<< ", " << ecls[nclu]<< ", " << extra_ << endl;
		
		
      thit[nclu]=1E6+1;
      for(Int_t ie=0;ie<ne[nclu];ie++)
	{
	  //get electron
	  //Double_t xel,yel,zel,tel,eel,dex,dey,dez;
		
		// Initial position [cm] and starting time [ns]
		double x0 = 1.1, y0 = 1.1, z0 = 0.;
		double t0 = 0.;
		// Initial energy [eV]
		double e0 = 100;
		// Initial direction
		// In case of a null vector, the initial direction is randomized.
		double dx0 = 0., dy0 = 0., dz0 = 0.;
		// Calculate an electron avalanche.
	    Tracking->GetElectron(ie,x0,y0,z0,t0,e0,dx0,dy0,dz0);
		//Tracking->GetElectron(ie,xel,yel,zel,tel,eel,dex,dey,dez);
	
		
	  //AvalMC drift
	  Double_t xw, yw, zw, tdrift;
	  aval->DriftElectron(x0,y0,z0,t0);
	  //aval->DriftElectron(xel,yel,zel,tel);
	  aval->GetDriftLinePoint(aval->GetNumberOfDriftLinePoints()-1, 
				  xw, yw, zw, tdrift);
	  // cout << "AvalMC (" << ie << ") " << xw<< ", " << yw << ", " << zw << ", " << tdrift << endl;
		//cout << "tdrift = " << tdrift<< endl;
	  if( (fabs(xw)>rCell/5.) ||  (fabs(yw)>rCell/5.) ){
	    //cout << "lost electron" << endl;
	    continue;
	  }
	  if(thit[nclu]>tdrift) thit[nclu] = tdrift;
	}
      if(thit[nclu]>1E6){ 
	//cout << "lost cluster" << endl; 
	continue;
      }
	
      dhit[nclu]=sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1]);
       //cout << "bdiff = " << (dhit[nclu]- sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1]))/dhit[nclu]<< endl;
       //cout << "dhit = " << dhit[nclu] << " =?= " << sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1])<< endl;
      nclu++;
      if(nclu>=MAXCLS) break;
		
		Int_t nion;
		aval->GetAvalancheSize(nelec[nclu], nion);
		cout << "electrons in aval: " << nelec[nclu] << " , ions in aval: " << nion << endl;
		ViewSignal* signalView = new ViewSignal();
		signalView->SetSensor(sensor);
		signalView->PlotSignal("s");
		
		TH1D* sig_hists = signalView->GetHistogram();
		integral[nclu] = sig_hists->Integral();
		
		cout << "signal integral: " << integral[nclu] << " fC = " << integral[nclu]/1.6E-4 << " electrons" << endl;
		gain[nclu]=integral[nclu]/1.6E-4;
		break;
    }

    tree->Fill();
  }


	
  sfile->cd();
  tree ->Write();
  sfile->Close();

  if(view) {
    ViewCell* view_cmp = new ViewCell();
    view_cmp->SetComponent(cmp); 
    view_cmp->Plot2d();
    view_drift->Plot();

    app->Run(true);
  }
}
