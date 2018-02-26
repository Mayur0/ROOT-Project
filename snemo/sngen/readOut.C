#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TDatime.h"
#include "TVector2.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRandom3.h"
#include <iostream>

using namespace std;

void readOut(){

	TChain *myTree = new TChain("clusters");
	TChain *myTree2 = new TChain("clusters");
	//myTree->AddFile("p=(e-)-E=1MeV-v=3000-rCell=2.20-he_90_Ib_10.gas-50clusters.root");
	myTree->AddFile("p=(mu-)-E=4000MeV-v=1250-rCell=2.20-he_90_Ib_10.gas-50clusters.root");
	myTree2->AddFile("p=(e-)-E=1MeV-v=2000-rCell=2.20-Ar_90_CO2_10.gas-100clusters.root");
	
	TGraph *g= new TGraph();
	TGraph *g1= new TGraph();
	TGraph *g2= new TGraph();
	TGraph *g3= new TGraph();
	TGraph *g4= new TGraph();
	TGraph *g5= new TGraph();
	TGraph *g6= new TGraph();
	TGraph *g7= new TGraph();
	TGraph *g8= new TGraph();
	TGraph *g9= new TGraph();
	TGraph *gt= new TGraph();

	g->SetNameTitle("g","xt relationship");
	g1->SetNameTitle("g1","xt relationship");
	g2->SetNameTitle("g2","xt relationship");
	g3->SetNameTitle("g3","xt relationship");
	g4->SetNameTitle("g4","xt relationship");
	g5->SetNameTitle("g5","xt relationship");
	g6->SetNameTitle("g6","xt relationship");
	g7->SetNameTitle("g7","xt relationship");
	g8->SetNameTitle("g8","xt relationship");
	g9->SetNameTitle("g9","xt relationship");
	gt->SetNameTitle("gt","minimum time distribution");

	Double_t dhit[200], thit[200], integral[200], gain[200];
	Int_t nclu, irow, icol, tmin;
	Int_t nelec[200];

	myTree->SetBranchAddress("nclu",&nclu);
	myTree->SetBranchAddress("irow",&irow);
	myTree->SetBranchAddress("icol",&icol);
	myTree->SetBranchAddress("dhit",&dhit[0]);
	myTree->SetBranchAddress("thit",thit);
	myTree->SetBranchAddress("integral",integral);
	myTree->SetBranchAddress("gain",gain);
	myTree->SetBranchAddress("nelec",&nelec);
	myTree->SetBranchAddress("tmin",&tmin);

	Int_t entries= myTree->GetEntries();

	for (int i=0; i<entries; i++)
	{
		myTree->GetEntry(i);
		for (int j=0; j<nclu; j++)
		{
			if (irow==0){
				g->SetPoint(g->GetN(),thit[j],dhit[j]);
			}
			else if (irow==1){
				g1->SetPoint(g1->GetN(),thit[j],dhit[j]);
			}
			else if (irow==2){
				g2->SetPoint(g2->GetN(),thit[j],dhit[j]);
			}
			else if (irow==3){
				g3->SetPoint(g3->GetN(),thit[j],dhit[j]);
			}
			else if (irow==4){
				g4->SetPoint(g4->GetN(),thit[j],dhit[j]);
			}
			else if (irow==5){
				g5->SetPoint(g5->GetN(),thit[j],dhit[j]);
			}
			else if (irow==6){
				g6->SetPoint(g6->GetN(),thit[j],dhit[j]);
			}
			else if (irow==7 && icol==0){
				g7->SetPoint(g7->GetN(),thit[j],dhit[j]);
			}
			else if (irow==8 && icol==1){
				g8->SetPoint(g8->GetN(),thit[j],dhit[j]);
			}
			else if (irow==7 && icol==1){
				g9->SetPoint(g9->GetN(),thit[j],dhit[j]);
			}
		}
	}


//	TGraph *ga= new TGraph();
//	ga->SetNameTitle("ga","xt relationship");
//	myTree2->SetBranchAddress("nclu",&nclu);
//	myTree2->SetBranchAddress("dhit",&dhit[0]);
//	myTree2->SetBranchAddress("thit",thit);
//	Int_t entries2= myTree2->GetEntries();
//	for (int i=0; i<entries2; i++)
//	{
//		myTree2->GetEntry(i);
//		for (int j=0; j<nclu; j++)
//		{
//			ga->SetPoint(ga->GetN(),thit[j],dhit[j]);
//		}
//	}

	TCanvas *c1 = new TCanvas("c1","Drift distances/times",200,10,700,500);

	c1->SetFillColor(42);
	c1->SetGrid();
	const Int_t m = nclu;

	//convert x from ns to microseconds by dividing by 1000
//	const Int_t n = 9000;
//	Double_t x[n], y[n];
//	for (Int_t i=0;i<n;i++) {
//		x[i] = i;
//		y[i] = (((x[i]/1000)*1.523)/((pow((x[i]/1000),0.572))+0.357));
//	}
//
//	TGraphErrors *gl = new TGraphErrors(n,x,y);

	g->SetMarkerColor(4);
	g1->SetMarkerColor(7);
	g2->SetMarkerColor(6);
	g3->SetMarkerColor(5);
	g4->SetMarkerColor(8);
	g5->SetMarkerColor(2);
	g6->SetMarkerColor(3);
	g7->SetMarkerColor(9);
	g8->SetMarkerColor(12);
	g9->SetMarkerColor(11);
	g->SetMarkerStyle(2);
	g1->SetMarkerStyle(2);
	g2->SetMarkerStyle(2);
	g3->SetMarkerStyle(2);
	g4->SetMarkerStyle(2);
	g5->SetMarkerStyle(2);
	g6->SetMarkerStyle(2);
	g7->SetMarkerStyle(2);
	g8->SetMarkerStyle(2);
	g9->SetMarkerStyle(2);
	g->SetTitle("Drift distances/times");
	g->GetXaxis()->SetLimits(0,6000);
	g->GetXaxis()->SetTitle("thit/ns");
	g->GetYaxis()->SetTitle("dhit/cm");
	g->Draw("AP");
	g1->Draw("P");
	g2->Draw("P");
	g3->Draw("P");
	g4->Draw("P");
	g5->Draw("P");
	g6->Draw("P");
	g7->Draw("P");
	g8->Draw("P");
	g9->Draw("P");

//	ga->SetMarkerColor(3);
//	ga->SetMarkerStyle(2);
//	ga->Draw("P");
//
//	gl->SetLineColor(2);
//	gl->SetLineWidth(2);

	TLegend *leg = new TLegend(0.6,0.1,0.9,0.5);
	leg->AddEntry("g","1st row","P");
	leg->AddEntry("g1","2nd row","p");
	leg->AddEntry("g2","3rd row","p");
	leg->AddEntry("g3","4th row","p");
	leg->AddEntry("g4","5th row","p");
	leg->AddEntry("g5","6th row","p");
	leg->AddEntry("g6","7th row","p");
	leg->AddEntry("g7","8th row","p");
	leg->AddEntry("g8","9th row","p");
	leg->AddEntry("g9","8th row, 2nd col","p");
	//leg->AddEntry("ga","electron","P");
	leg->Draw();
	
	//gl->Draw("L");
	
	//TH1F *hist = new TH1F("hist","Gain integral",10000,0,10000);
	//for (Int_t indx=0;indx<myTree->GetEntries();indx++){
		//myTree->GetEntry(indx);
		//for (int j=0; j<nclu; j++)
		//{
			//hist->Fill(nelec[j]);
		//}
	//}
//	hist->GetXaxis()->SetTitle("Gain");
//	hist->GetYaxis()->SetTitle("Frequency");
//	hist->Fit("gaus");
//	hist->Draw();

//	TProfile *hprof = new TProfile("hprof","Drift histogram",40,0,6000,0,3.5);
//	//bins,xmin,xmax,ymin,ymax
//	for (Int_t indx=0;indx<myTree->GetEntries();indx++){
//		myTree->GetEntry(indx);
//		for (int j=0; j<nclu; j++)
//		{
//			hprof->Fill(thit[j],dhit[j]);
//		}
//	}
//
//	hprof->GetXaxis()->SetTitle("thit/ns");
//	hprof->GetYaxis()->SetTitle("dhit/cm");
//	TF1 *func = new TF1("func","(((x/1000)*[0])/((pow((x/1000),[1]))+[2]))",0,5000);
//	func->SetParameters(2.3997,0.7731,0.3347);
//	hprof->Fit("func","","",0,2600);
//
//	hprof->Draw();
//
//	TLegend *hleg = new TLegend(0.6,0.1,0.9,0.4);
//	hleg->AddEntry("","A = 1.532","");
//	hleg->AddEntry("","B = 0.465","");
//	hleg->AddEntry("","C = 0.228","");
//	hleg->Draw();
	
	Double_t tmin_0, tmin_1, tmin_8, dmin_0, dmin_1, dmin_8, smin;

	TH1F *hist_dmin= new TH1F("hist_dmin","minimum distance",50,-2,2);

	for (int i=0; i<entries; i++){
		myTree->GetEntry(i);
		
//		dmin[] = (((tmin[]/1000)*1.532)/((pow((tmin[]/1000),0.465))+0.228));

		if (irow==0 && icol==0){
			tmin_0 = tmin;
			cout << "tmin_0 = " << tmin_0 << endl;
			dmin_0 = (((tmin_0/1000.)*1.532)/((pow((tmin_0/1000.),0.465))+0.228));
			cout << "dmin_0 = " << dmin_0 << endl;
		}
		else if (irow==1 && icol==0){
			tmin_1 = tmin;
			cout << "tmin_1 = " << tmin_1 << endl;
			dmin_1 = (((tmin_1/1000.)*1.532)/((pow((tmin_1/1000.),0.465))+0.228));
			cout << "dmin_1 = " << dmin_1 << endl;
		}
		else if (irow==8 && icol==1){
			tmin_8 = tmin;
			cout << "tmin_8 = " << tmin_8 << endl;
			dmin_8 = (((tmin_8/1000.)*1.532)/((pow((tmin_8/1000.),0.465))+0.228));
			cout << "dmin_8 = " << dmin_8 << endl;
			
			smin = ((dmin_8+dmin_0)/2)-dmin_1;
			if(smin==0){
				continue;
			}
			cout << "smin = " << smin << endl;
			hist_dmin->Fill(smin);

		}
		
		
	}
	hist_dmin->GetXaxis()->SetTitle("dmin/cm");
	hist_dmin->GetYaxis()->SetTitle("Frequency");
	hist_dmin->Draw();

	
}
