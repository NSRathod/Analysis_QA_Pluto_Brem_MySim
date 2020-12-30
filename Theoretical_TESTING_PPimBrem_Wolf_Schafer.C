#include "TPaveText.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TText.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TMatrixD.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <vector>
//#include <random>
#include "TCanvas.h"
#include "TPad.h"
#include "TCutG.h"
#include <cmath>
#include <TFile.h>
#include <TLorentzVector.h>

#define SQR(x) ((x)*(x))
using namespace std;
void normalize(TH1* hist);
double openingangle(const TLorentzVector& a, const TLorentzVector& b)
{
  return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Vect().Mag() * b.Vect().Mag() ) );
}

void Theoretical_TESTING_PPimBrem_Wolf_Schafer()
{
ifstream f1a;
//f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/PimP_Delta_dalitz_P690_1M.evt");
//f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PimDelta_P690MeV.evt");
//f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PP_NEW_Bremsstrahlung.evt");
f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PPim_Bremsstrahlung_P690MeV.evt");
//f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_PPim_brems_new/p_ee_0.evt");

 const float PI= TMath::Pi();

 double pion_momentum = 0.690;//GeV/c
 const float mpi=0.13957;

 double pion_energy = sqrt(pion_momentum*pion_momentum + mpi*mpi);

 const float mp=0.93827;
 const float me=0.510998*0.001;

 const float alpha = 1./137.;
 const float sqrts=1.500;
 const float sigma_el=18.13;//mb (SAID)
 //const float sigma_el=18.13*1000.;//microbarn (SAID)
 //const float A=alpha*alpha/(6.*PI*PI*PI);
 const float A=alpha*alpha*1.9/(3.*PI*PI*PI);
 int flag1=0, flag2=0;

 
  TLorentzVector *proj;
  proj = new TLorentzVector(0,0,pion_momentum, pion_energy); //---------------PION BEAM
  TLorentzVector *targ;
  targ = new TLorentzVector(0,0,0, 0.93827231);          //-------------------PROTON TARGET
  TLorentzVector *beam;
  beam = new TLorentzVector(0,0,0,0);
  *beam = *proj + *targ;

  TLorentzVector *p;
  p = new TLorentzVector(0,0,0,0);            //---------------PION BEAM
  TLorentzVector *pim;
  pim = new TLorentzVector(0,0,0,0);          //-------------------PROTON TARGET
  TLorentzVector *ep;
  ep = new TLorentzVector(0,0,0,0);
  TLorentzVector *em;
  em = new TLorentzVector(0,0,0,0);
  TLorentzVector *ep1;
  ep1 = new TLorentzVector(0,0,0,0);
  TLorentzVector *em1;
  em1 = new TLorentzVector(0,0,0,0);
  TLorentzVector *mm;
  mm = new TLorentzVector(0,0,0,0);
  TLorentzVector *invm;
  invm = new TLorentzVector(0,0,0,0);
  TLorentzVector *InvCM;
  InvCM = new TLorentzVector(0,0,0,0);
		   
  Double_t gloW; 

   TH1F *epemInvM = new TH1F("epemInvM","epemInvM",100,0.,1.);
   epemInvM->Sumw2();
   TH1F *epemInvMcal = new TH1F("epemInvMcal","epemInvMcal",100,0.,1.);
   epemInvMcal->Sumw2();
   TH1F *hepemMM = new TH1F("hepemMM","hepemMM",100,0.,1.2);
   hepemMM->Sumw2();
   TH1F *hOA = new TH1F("hOA","hOA",180, 0.0, 180);
   hOA->Sumw2();
   TH1F *hydil = new TH1F("hydil","hydil",100, -2.5, 2.5);
   hydil->Sumw2();
   TH1F *hptdil = new TH1F("hptdil","hptdil",100, 0, 0.5);
   hptdil->Sumw2();
   
   TH2F *sigInvM = new TH2F("sigInvM","sigInvM",100,0.,0.4,100,0,3);
   TH2F *yInvM = new TH2F("yInvM","yInvM",100,0.,0.4,100,-2.5,2.5);
   TH2F *ptInvM = new TH2F("ptInvM","ptInvM",100,0.,0.4,100,0,0.5);

   TH1F* epem_4pi_Wolf_Schafer = new TH1F("epem_4pi_Wolf_Schafer","epem_4pi_Wolf_Schafer",111,0.0,0.555);
   epem_4pi_Wolf_Schafer->Sumw2();

   TH1F* epem_4pi_phasespace = new TH1F("epem_4pi_phasespace","epem_4pi_phasespace",111,0.0,0.555);
   epem_4pi_phasespace->Sumw2();
   
   
   float sig[8]={0.};
   float nb[8]={0.};
   float nbEv=1000000.;


   ofstream f1;
   f1.open ("PPim_Bremsstrahlung_evt_wts.evt");
   
while (!f1a.eof()){


  float p1a,wa,w1a,w2a,w3a;
  float a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13;
  float b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13;
  float c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13;
  float d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13;

  f1a>>p1a>>wa>>w1a>>w2a>>w3a;
  f1a>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12>>a13;
  f1a>>b1>>b2>>b3>>b4>>b5>>b6>>b7>>b8>>b9>>b10>>b11>>b12>>b13;
  f1a>>c1>>c2>>c3>>c4>>c5>>c6>>c7>>c8>>c9>>c10>>c11>>c12>>c13;
  f1a>>d1>>d2>>d3>>d4>>d5>>d6>>d7>>d8>>d9>>d10>>d11>>d12>>d13;

  // cout<<a1<<endl;
	 
  TVector3 m3,vert, p1, p2, p3, p4,p5,p6;
  TLorentzVector *lv = new TLorentzVector;
  TLorentzVector *Theta_CM;
  TLorentzVector *lv1 = new TLorentzVector;
  TLorentzVector *lv2 = new TLorentzVector;	
  TLorentzVector *lv3 = new TLorentzVector;
  TLorentzVector *lv4 = new TLorentzVector;
  
  lv->SetPxPyPzE(c2,c3,c4,c1);	 
  lv1->SetPxPyPzE(d2,d3,d4,d1);	 

  lv2->SetPxPyPzE(c2,c3,c4,c1);
  lv2->Boost(-(*beam).BoostVector());
  lv3->SetPxPyPzE(d2,d3,d4,d1);
  lv3->Boost(-(*beam).BoostVector());

  flag1=0;
  flag2=0;
  
    if(c9==2/* && th>18 && th<85*/){
    p3.SetXYZ(lv->Px(),lv->Py(),lv->Pz());
    ep->SetVectM(p3, me);
    flag1=1;//
    }		      
  
    if(d9==3/* && th>18 && th<85*/){
    p4.SetXYZ(lv1->Px(),lv1->Py(),lv1->Pz());
    em->SetVectM(p4, me);
    flag2=1;	      
    }  
  
    if(c9==2/* && th>18 && th<85*/){
    p5.SetXYZ(lv2->Px(),lv2->Py(),lv2->Pz());
    ep1->SetVectM(p5, me);
    flag1=1;//
    }		      
  
    if(d9==3/* && th>18 && th<85*/){
    p6.SetXYZ(lv3->Px(),lv3->Py(),lv3->Pz());
    em1->SetVectM(p6, me);
    flag2=1;	      
    }  

 
 
if(flag1 && flag2){
   *mm=*beam-*ep-*em;
   *invm = *ep + *em;
   *InvCM = *ep1 + *em1;
   float Invmass = invm->M();
   float InvEng = invm->E();
   float InvMom = invm->P();
   float pt = invm->Pt();
   float y = invm->Rapidity();
   float y_ep = ep->Rapidity();
   float y_em = em->Rapidity();

   float pt_ep = ep->Pt();
   float pt_em = em->Pt();
   float missmass = mm->M();

   float InvEngCM = InvCM->E(); 

   //cout<<InvEng<<"    "<<InvEngCM<<endl;
   //cout<<"pt_em: "<<pt_em<<" pt_ep: "<<pt_ep<<" pt_dil: "<<pt<<" y_dil: "<<y<<" y_ep: "<<y_ep<<" y_em: "<<y_em<<endl;
   //cout<<" pt_dil: "<<pt<<" y_dil: "<<y<<endl;
   
   //cout<<event<<"\t\t\t"<<"Mass: " <<Invmass<<"\t\t\t"<<"Energy: "<<InvEng<<"\t\t\t"<<"Momentum: "<<InvMom<<endl;
   Double_t Opang = openingangle(*ep,*em);
   Double_t Opang_deg = Opang*TMath::RadToDeg();

   //cout<<Invmass<<endl;
   
   /*   float s2=sqrts*sqrts+Invmass*Invmass-(2.*InvEngCM*sqrts);
   float R2s=sqrt(1.-((mpi+mp)*(mpi+mp))/(sqrts*sqrts));
   float R2s2=sqrt(1.-(mpi+mp)*(mpi+mp)/s2);

   float sigma_quasiel= ((sqrts*sqrts-(mpi+mp)*(mpi+mp))/(2.*mpi*mpi))*sigma_el;

   float dsigma=A*(sigma_quasiel/(Invmass*InvEng*InvEng))*(R2s2/R2s)*y*pt*pt;*/
   float AA = 1;
   //float s2=(sqrts*sqrts)+(Invmass*Invmass)-(2.*InvEngCM*sqrts);
   float s2=(sqrts*sqrts)+(Invmass*Invmass)-(2.*Invmass*sqrts);
   float R2s=sqrt(1.-(((mpi+mp)*(mpi+mp))/(sqrts*sqrts)));
   float R2s2=sqrt(1.-(((mpi+mp)*(mpi+mp))/s2));
   float sigma_quasiel= ((sqrts*sqrts-((mpi+mp)*(mpi+mp)))/(2.*mpi*mpi))*sigma_el;
   //float sigma_quasiel= ((sqrts*sqrts/(4*(mpi*mpi)))-1)*2*sigma_el;

   //float dsigma=A*(sigma_quasiel/(Invmass*InvEng*InvEng))*(R2s2/R2s)*y*pt*pt*(1/(Invmass));    //standard approach
   float dsigma=A*(sigma_quasiel/(Invmass))*(R2s2/R2s)*log((sqrts-(mpi+mp))/Invmass);            //standard approach
   //cout<<Invmass<<" "<<dsigma<<endl;
   //if (TMath::IsNaN(dsigma))continue;
   if((((mpi+mp)*(mpi+mp))/(s2)) > 1)continue;
   //if(dsigma<0.01)
   //cout<<Invmass<<" "<<dsigma<<endl;
   //hOA->Fill(Opang_deg);
   epemInvMcal->Fill(Invmass,dsigma);
   epemInvM->Fill(Invmass,dsigma*(1./nbEv));
   epem_4pi_Wolf_Schafer->Fill(Invmass,dsigma*(1./nbEv));
   epem_4pi_phasespace->Fill(Invmass,(1./nbEv));
   //epemInvM->Fill(Invmass);
   sigInvM->Fill(Invmass,dsigma);
   
   yInvM->Fill(Invmass,y);
   ptInvM->Fill(Invmass,pt);
 
   hydil->Fill(y);
   hptdil->Fill(pt);

   
   //if(y>0 && y<0.5 && pt>0.1 && pt<0.2){
     if (Invmass>0 && Invmass<=0.05){ sig[0]+=dsigma;nb[0]++;}
     else if (Invmass>0.05 && Invmass<=0.1){sig[1]+=dsigma;nb[1]++;}
     else if (Invmass>0.1 && Invmass<=0.15){sig[2]+=dsigma;nb[2]++;}
     else if (Invmass>0.15 && Invmass<=0.2){sig[3]+=dsigma;nb[3]++;}
     else if (Invmass>0.2 && Invmass<=0.25){sig[4]+=dsigma;nb[4]++;}
     else if (Invmass>0.25 && Invmass<=0.3){sig[5]+=dsigma;nb[5]++;}
     else if (Invmass>0.3 && Invmass<=0.35){sig[6]+=dsigma;nb[6]++;}
     else if (Invmass>0.35 && Invmass<=0.4){sig[7]+=dsigma;nb[7]++;}

//////////////////////////////////////////////////////////////////////////////////////////////////////
/*                                                                                                                                                                                  
cout<<p1a<<" "<<wa<<" "<<w1a<<" "<<w2a<<" "<<w3a<<endl;                                                                                                                              
cout<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a13*dsigma<<endl;                                        
cout<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<" "<<b11<<" "<<b12<<" "<<b13*dsigma<<endl;                                        
cout<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<" "<<c11<<" "<<c12<<" "<<c13*dsigma<<endl;                                        
cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<" "<<d7<<" "<<d8<<" "<<d9<<" "<<d10<<" "<<d11<<" "<<d12<<" "<<d13*dsigma<<endl;                                        
*/

  f1<<p1a<<" "<<wa<<" "<<w1a<<" "<<w2a<<" "<<w3a<<endl;
  f1<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<dsigma<<"\n";
  f1<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<" "<<b11<<" "<<b12<<" "<<dsigma<<"\n";
  f1<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<" "<<c11<<" "<<c12<<" "<<dsigma<<"\n";
  f1<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<" "<<d7<<" "<<d8<<" "<<d9<<" "<<d10<<" "<<d11<<" "<<d12<<" "<<dsigma<<"\n";   
//////////////////////////////////////////////////////////////////////////////////////////////////////
     }
 }

 f1.close();
 
 float sigfin[8]={0.};
 float invM[8]={0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375};
 float dM=0.05;
 float dy=0.5;
 float dpt=0.1;
 
 for (int i=0;i<8;i++){

   //sigfin[i]=sig[i]/(nb[i]*dM*dy*dpt);
   sigfin[i]=sig[i]/(nb[i]*dM);
   cout<<i<<"  "<<nb[i]<<"   "<<invM[i]<<"  "<<sigfin[i]<<endl;

 }

 epemInvMcal->Scale(1./1000000.);
 
   normalize(epemInvMcal);
   normalize(epemInvM);
   normalize(epem_4pi_Wolf_Schafer);
   normalize(epem_4pi_phasespace);
 epemInvMcal->Draw("same,hist");

 
 TGraph *gr=new TGraph(8,invM,sigfin);
 gr->SetMarkerStyle(21);
 gr->SetMarkerSize(1.);
 
 TCanvas *can1=new TCanvas("brems_theory","brems_theory");
 can1->cd ();
 gPad->SetLogy();
 gr->Draw("al");
 
 
 
 TCanvas *can5=new TCanvas("brems_hist","brems_hist");
 //can->Divide(3,2);
 

 can5->cd();
 gPad->SetLogy();
 epem_4pi_Wolf_Schafer->SetLineWidth(2);
 epem_4pi_Wolf_Schafer->SetLineColor(kRed);
 epem_4pi_phasespace->SetLineWidth(2);
 epem_4pi_phasespace->SetLineColor(kBlue);
 //epemInvMcal->Draw("same,hist");
 epem_4pi_Wolf_Schafer->Draw("same");
 epem_4pi_phasespace->Draw("same"); 
 //epemInvM->SetLineWidth(2);
 //epemInvM->Draw("same,hist");
 /* 
 TCanvas *can=new TCanvas("brems_hist","brems_hist",1000,600);
 can->Divide(2,1);
 can->cd(1);
 gPad->SetLogy();
 normalize(epemInvM);
 epemInvM->Draw("same");

 can->cd(2);
 gPad->SetLogy();
 gr->Draw("alsame");
 gr2->Draw("psame");
*/ 

 TFile *output = new TFile("Wolf_Schafer_ph_Inv_PimP_1GeV_Bremsstrahlung_CHECK.root","RECREATE");
 output->cd();
 //CSvsInvM->Write();
 //can->Write();
 can1->Write();
 can5->Write();
 epemInvMcal->Write();
 epemInvM->Write();
 epem_4pi_Wolf_Schafer->Write();
 epem_4pi_phasespace->Write();
 //can6->Write();
 output->Write();
 output->Close();
 
 /*
 can->cd(2);
 sigInvM->Scale(1./sigInvM->Integral());
 
 sigInvM->Draw("colz"); 

 can->cd(3);
 hydil->Draw();
 
 can->cd(4);
 hptdil->Draw();


 can->cd(5);
 yInvM->Draw("colz"); 

 can->cd(6);
 ptInvM->Draw("colz"); 
 */
 
} 
 
void normalize(TH1* hist)
{
  for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
    {
      hist->SetBinContent( j, hist->GetBinContent(j) / hist->GetBinWidth(j) );
      //         hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
      hist->SetBinError( j, hist->GetBinError(j) / hist->GetBinWidth(j) );
    }
}




