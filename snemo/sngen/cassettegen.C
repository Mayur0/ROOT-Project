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

#define DEBUG false

using namespace Garfield;
using namespace std;

int main(int argc,char *argv[]){

  //TString gasFile = "../data/gas/he_95_ar_1_eth_4.gas";
  TString gasFile = "../data/gas/He_95_Ar_1_Eth_4.gas";
  int Ntrk=2;
  Bool_t view = false;
  // Cell radius [cm]
  float rCell = 2.2;
  float WireGap =  0.9115;
  const float XWallGap =  1.4;
  const float CaloGap =  1;
  // Wire radius [cm]
  const double rSWire = 5.e-3; //50 um *FIXME* this is actually a diameter
  const double rFWire = 5.e-3; //50 um *FIXME* this is actually a diameter
  // Half-length of the tube [cm]
  const double lTube = 360.;
  // Voltages
  const double vWire = 1250.;

  const int ncols =1;
  const int nrows =2;

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
       WireGap = f*WireGap;
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

  float Xend = XWallGap+(ncols*rCell*2);
  float Yend = CaloGap+ (nrows*rCell*2);

  //Make a large tube of gas including all the chamber
  SolidTube *tube = new SolidTube(Xend/2.,Yend/2.,0.,rSWire,Xend*sqrt2,lTube);
  geo->AddSolid(tube,gas);

  cmp->SetGeometry(geo);
  //Add ground planes
  //AddPlaneX (const double x, const double voltage, const std::string label)
  cmp->AddPlaneX (0, 0, "x wall");
  cmp->AddPlaneY (0, 0, "calo wall");

  cmp->AddPlaneX (Xend+XWallGap, 0, "x wall2");
  cmp->AddPlaneY (Yend+CaloGap, 0, "calo wall2");

  //cmp->SetPeriodicityX(ncols*rCell*2);
  //cmp->SetPeriodicityY(nrows*rCell*2);

  //new cohordinate system, centered in the calo corner
  for(int col=0; col<ncols; col++){
    float Xorig = XWallGap+(col*rCell*2);
    for(int row=0; row<nrows; row++){
      float Yorig = CaloGap+(row*rCell*2);

      //Add sense wires
      cmp->AddWire(Xorig+rCell,Yorig+rCell,rSWire,vWire,Form("s%d%d",col,row)); 
      //Add field wires
      cmp->AddWire(Xorig,Yorig+rCell,rFWire,0.,"f");
      cmp->AddWire(Xorig,Yorig+rCell-WireGap,rFWire,0.,"f");
      cmp->AddWire(Xorig,Yorig+rCell+WireGap,rFWire,0.,"f");

      cmp->AddWire(Xorig+rCell,Yorig,rFWire,0.,"f");
      cmp->AddWire(Xorig+rCell-WireGap,Yorig,rFWire,0.,"f");
      cmp->AddWire(Xorig+rCell+WireGap,Yorig,rFWire,0.,"f");

    }

    cmp->AddWire(Xorig+rCell,Yend,rFWire,0.,"f");
    cmp->AddWire(Xorig+rCell-WireGap,Yend,rFWire,0.,"f");
    cmp->AddWire(Xorig+rCell+WireGap,Yend,rFWire,0.,"f");
  }


  for(int row=0; row<nrows; row++){
    float Yorig = CaloGap+(row*rCell*2);
    cmp->AddWire(Xend,Yorig+rCell,rFWire,0.,"f");
    cmp->AddWire(Xend,Yorig+rCell-WireGap,rFWire,0.,"f");
    cmp->AddWire(Xend,Yorig+rCell+WireGap,rFWire,0.,"f");    
  }

#define  MAXCLS 10
	
	Int_t nclu, ne[MAXCLS];
	Double_t cls[MAXCLS][3], tcls[MAXCLS], ecls[MAXCLS], dhit[MAXCLS], thit[MAXCLS], dhit_2[MAXCLS], thit_2[MAXCLS];
	Double_t orig[3],  dir[3], p1[3], p2[3];
	Double_t b, theta, phi;
	
	TFile *sfile = new TFile(title+Form("-%dclusters.root",Ntrk),"RECREATE");
	cout << sfile->GetName() << endl;
	TTree *tree = new TTree("clusters","clusters");
	tree->Branch("b",&b,"b/D");
	tree->Branch("theta",&theta,"theta/D");
	tree->Branch("phi",&phi,"phi/D");
	tree->Branch("nclu",&nclu,"nclu/I");
	
	tree->Branch("ne",ne,"ne[nclu]/I");
	tree->Branch("cls",cls,"cls[nclu][3]/D");
	tree->Branch("tcls",tcls,"tcls[nclu]/D");
	tree->Branch("ecls",ecls,"ecls[nclu]/D");
	tree->Branch("dhit",dhit,"dhit[nclu]/D");
	tree->Branch("thit",thit,"thit[nclu]/D");
	tree->Branch("dhit_2",dhit_2,"dhit_2[nclu]/D");
	tree->Branch("thit_2",thit_2,"thit_2[nclu]/D");
	
	tree->Branch("orig",orig,"orig[3]/D");
	tree->Branch("dir",dir,"dir[3]/D");
	tree->Branch("p1",p1,"p1[3]/D");
	tree->Branch("p2",p2,"p2[3]/D");
	
	Sensor* sensor = new Sensor();
	AvalancheMC* aval = new AvalancheMC();
	TrackHeed *Tracking = new TrackHeed();
	
//	TString sname= Form("s%d%d",0,0);
//	cmp->AddReadout(sname.Data());
//	sname= Form("s%d%d",0,1);
//	cmp->AddReadout(sname.Data());
	
	for(Int_t itrk=0;itrk<Ntrk;itrk++){
		// cout << endl << "event " << itrk << endl;
//		if(itrk%10==0)
		{ cout << "track n. = " << itrk << endl;}
//		TVector3 vorig, vdir;
//		TVector3 vp1,vp2,vzero(0,0,0);
//
//		b = gRandom->Uniform(0,sqrt2*rCell);
//		phi=gRandom->Uniform(0,TMath::TwoPi());
//		theta=TMath::PiOver2();
//
//		getTrack2D(b,phi,rCell,vp1,vp2);
//		if( (vp1==vzero) | (vp2==vzero)){ itrk--; continue;}
//
//		vorig=vp1;
//		vdir=vp2-vp1;
//		vdir.SetTheta(theta);
//
//		vorig.GetXYZ(orig);
//		vdir.GetXYZ(dir);
//		vp1.GetXYZ(p1);
//		vp2.GetXYZ(p2);
		
		//    Tracking->NewTrack(vorig.X(),vorig.Y(),vorig.Z(), //Position
		//		       0.,                   //Time
		//		       vdir.X(),vdir.Y(),vdir.Z());  //Direction
		
  for(int irow=0; irow<nrows; irow++){
	  for(int icol=0; icol<ncols; icol++){
	  cout << "row = " << irow << endl;
	  cout << "col = " << icol << endl;
	  
	  TString sname= Form("s%d%d",icol,irow);
	  cmp->AddReadout(sname.Data());
	  
	  Double_t Xorig = XWallGap+(icol*rCell*2);
	  Double_t Yorig = CaloGap+(irow*rCell*2);
	  Double_t wire_x = Xorig+rCell, wire_y = Yorig+rCell;

	  //Define a sensor
	  
	  sensor->AddComponent(cmp);
	  
	  sensor->AddElectrode(cmp,sname.Data());
	  sensor->SetTimeWindow(0., 0.1, 1000);
	  sensor->SetArea(0,0,-lTube,ncols*2*rCell,nrows*2*rCell,lTube);
	  
	  //Tracking and avalanches

	  aval->SetSensor(sensor);
	  aval->EnableDiffusion();
  
	  Tracking->SetParticle(particle.Data());
	  Tracking->SetEnergy(energy);
	  Tracking->SetSensor(sensor);
	  
//	  ViewDrift* view_drift = 0;
//	  if(view){
//		  view_drift = new ViewDrift();
//		  view_drift->SetArea(0,0,-lTube,ncols*2*rCell,nrows*2*rCell,lTube);
//		  Tracking->EnablePlotting(view_drift);
//	  }
  
	  Tracking->NewTrack(XWallGap+0.5*rCell,0,0, //Position
						 0.,                   //Time
						 0,1,0);  //Direction

		Double_t extra_;
	
    nclu = 0;
	  
	if(irow==0){
    while(Tracking->GetCluster(cls[nclu][0],cls[nclu][1],cls[nclu][2],
			       tcls[nclu],
			       ne[nclu],ecls[nclu],extra_)){
		cout << nclu << " - " << cls[nclu][0] << ", " << cls[nclu][1]<< ", " << cls[nclu][2]<< ", " << tcls[nclu]<< ", " << ne[nclu]<< ", " << ecls[nclu]<< ", " << extra_ << endl;
		
			
      thit[nclu]=1E6+1;
      for(Int_t ie=0;ie<ne[nclu];ie++)
	{
	  //get electron
	  Double_t xel,yel,zel,tel,eel,dex,dey,dez;
	  Tracking->GetElectron(ie,xel,yel,zel,tel,eel,dex,dey,dez);
   
	  //AvalMC drift
	  Double_t xw, yw, zw, tdrift;
	  aval->DriftElectron(xel,yel,zel,tel);
	  aval->GetDriftLinePoint(aval->GetNumberOfDriftLinePoints()-1, 
				  xw, yw, zw, tdrift);
	  // cout << "AvalMC (" << ie << ") " << xw<< ", " << yw << ", " << zw << ", " << tdrift << endl;
		
		if( (fabs(xw-wire_x)>rCell/5.) ||  (fabs(yw-wire_y)>rCell/5.) ){
			if(DEBUG){cout << "lost electron" << endl;}
			continue;
		}
		
		if(thit[nclu]>tdrift) thit[nclu] = tdrift;
	}
		
		if(thit[nclu]>1E6){
	cout << "lost cluster" << endl;
	continue;
      } 

      dhit[nclu]=sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1]);
       //cout << "bdiff = " << (dhit[nclu]- sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1]))/dhit[nclu]<< endl;
		cout << "dhit = " << dhit[nclu] << " =?= " << sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1])<< endl;
      nclu++;
      if(nclu>=MAXCLS) break;
    	}
	  
	  }
	else if(irow==1){
		  while(Tracking->GetCluster(cls[nclu][0],cls[nclu][1],cls[nclu][2],
									 tcls[nclu],
									 ne[nclu],ecls[nclu],extra_)){
			  cout << nclu << " - " << cls[nclu][0] << ", " << cls[nclu][1]<< ", " << cls[nclu][2]<< ", " << tcls[nclu]<< ", " << ne[nclu]<< ", " << ecls[nclu]<< ", " << extra_ << endl;
			  
			  
			  thit_2[nclu]=1E6+1;
			  for(Int_t ie=0;ie<ne[nclu];ie++)
			  {
				  //get electron
				  Double_t xel,yel,zel,tel,eel,dex,dey,dez;
				  Tracking->GetElectron(ie,xel,yel,zel,tel,eel,dex,dey,dez);
				  
				  //AvalMC drift
				  Double_t xw, yw, zw, tdrift;
				  aval->DriftElectron(xel,yel,zel,tel);
				  aval->GetDriftLinePoint(aval->GetNumberOfDriftLinePoints()-1,
										  xw, yw, zw, tdrift);
				  // cout << "AvalMC (" << ie << ") " << xw<< ", " << yw << ", " << zw << ", " << tdrift << endl;
				  
				  if( (fabs(xw-wire_x)>rCell/5.) ||  (fabs(yw-wire_y)>rCell/5.) ){
					  if(DEBUG){cout << "lost electron" << endl;}
					  continue;
				  }
				  
				  if(thit_2[nclu]>tdrift) thit_2[nclu] = tdrift;
			  }
			  
			  if(thit_2[nclu]>1E6){
				  cout << "lost cluster" << endl;
				  continue;
			  }
			  
			  dhit_2[nclu]=sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1]);
			  //cout << "bdiff = " << (dhit_2[nclu]- sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1]))/dhit_2[nclu]<< endl;
			  cout << "dhit = " << dhit_2[nclu] << " =?= " << sqrt(cls[nclu][0]*cls[nclu][0]+cls[nclu][1]*cls[nclu][1])<< endl;
			  nclu++;
			  if(nclu>=MAXCLS) break;
		  }
		  
	  }
	
    tree->Fill();
	  }
  }
  }
  sfile->cd();
  tree ->Write();
  sfile->Close();

  if(view) {
    ViewCell* view_cmp = new ViewCell();
    view_cmp->SetComponent(cmp); 
    view_cmp->Plot2d();

    //view_drift->Plot();
    app->Run(true);
  }
}
