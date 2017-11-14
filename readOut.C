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
	//TChain *myTree2 = new TChain("clusters");
	//myTree->AddFile("p=(e-)-E=1MeV-v=1830-rCell=2.20-He_95_Ar_1_Eth_4.gas-100clusters.root");
	myTree->AddFile("p=(e-)-E=1MeV-v=1830-rCell=2.20-Ar_90_CO2_10.gas-100clusters.root");
	
	TGraph *g= new TGraph();
	g->SetNameTitle("g","xt relationship");
	
	Double_t dhit[50], thit[50];
	Int_t nclu;
	
	myTree->SetBranchAddress("nclu",&nclu);
	myTree->SetBranchAddress("dhit",&dhit[0]);
	myTree->SetBranchAddress("thit",thit);
	Int_t entries= myTree->GetEntries();
	for (int i=0; i<entries; i++)
	{
		myTree->GetEntry(i);
		for (int j=0; j<nclu; j++)
		{
			g->SetPoint(g->GetN(),thit[j],dhit[j]);
		}
	}
	
	//TGraph *g2= new TGraph();
	//g2->SetNameTitle("g2","xt relationship");
	//myTree2->SetBranchAddress("nclu",&nclu);
	//myTree2->SetBranchAddress("dhit",&dhit[0]);
	//myTree2->SetBranchAddress("thit",thit);
	//Int_t entries2= myTree2->GetEntries();
	//for (int i=0; i<entries2; i++)
	//{
		//myTree2->GetEntry(i);
		//for (int j=0; j<nclu; j++)
		//{
			//g2->SetPoint(g2->GetN(),thit[j],dhit[j]);
		//}
	//}
	
	TCanvas *c1 = new TCanvas("c1","Drift distances/times",200,10,700,500);
	
	c1->SetFillColor(42);
	c1->SetGrid();
	const Int_t m = nclu;
	
	const Int_t n = 9000;
	Double_t x[n], y[n];
	for (Int_t i=0;i<n;i++) {
		x[i] = i;
		y[i] = (((x[i]/1000)*1.523)/((pow((x[i]/1000),0.572))+0.357));
	}
	
	TGraphErrors *gl = new TGraphErrors(n,x,y);
	
	g->SetMarkerColor(4);
	g->SetMarkerStyle(2);
	g->SetTitle("Drift distances/times");
	g->GetXaxis()->SetTitle("thit/ns");
	g->GetYaxis()->SetTitle("dhit/cm");
	g->Draw("AP");
	//g2->SetMarkerColor(4);
	//g2->SetMarkerStyle(2);
	//g2->Draw("P");
	
	gl->SetLineColor(2);
	gl->SetLineWidth(2);
	
	//TLegend *leg = new TLegend(0.1,0.8,0.48,0.9);
	//leg->AddEntry("g","Ar_90_CO2_10","P");
	//leg->AddEntry("g2","He_95_Ar_1_Eth_4","P");
	//leg->Draw();
	
	gl->Draw("L");

	TProfile *hprof = new TProfile("hprof","Drift histogram",100,0,10000,0,3.5);
	//100 bins,xmin,xmax,ymin,ymax
	for (Int_t indx=0;indx<myTree->GetEntries();indx++){
		myTree->GetEntry(indx);
		for (int j=0; j<nclu; j++)
		{
			hprof->Fill(thit[j],dhit[j]);
		}
	}
	
	hprof->GetXaxis()->SetTitle("thit/ns");
	hprof->GetYaxis()->SetTitle("dhit/cm");
	TF1 *func = new TF1("func","(((x/1000)*[0])/((pow((x/1000),[1]))+[2]))",0,3000);
	func->SetParameters(1.604,0.655,0.437);
	hprof->Fit("func","","",0,2600);
	//hprof->Fit("pol5")
	hprof->Draw();
	
}
