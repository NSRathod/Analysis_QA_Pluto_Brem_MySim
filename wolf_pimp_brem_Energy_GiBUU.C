#include "TPaveText.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
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
#include <TLorentzVector.h>
#include <TFile.h>


#define SQR(x) ((x)*(x))
using namespace std;

double openingangle(const TLorentzVector& a, const TLorentzVector& b)
{
  return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Vect().Mag() * b.Vect().Mag() ) );
}
void normalize(TH1* hist);
void wolf_pimp_brem_Energy_GiBUU(){
  ifstream f1a ;
  f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PPim_Bremsstrahlung_P690MeV.evt");
  
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_PPim_brems_new/p_ee_0.evt");

 const float PI= TMath::Pi();

 double p_momentum =  0.690;//GeV/c
 const float mpi=0.13957;
 const float mp=0.93827;
 const float me=0.510998*0.001;

 double p_energy = sqrt(p_momentum*p_momentum + mpi*mpi);

 const float alpha = 1./137.;
 const float sqrts=1.5;
 const float sigma_el=18.30;//mb 
 const float A=alpha*alpha/(6.*PI*PI*PI);
 const float A1=2*alpha*alpha/(3.*PI*PI);
 int flag1=0, flag2=0;
 
  TLorentzVector *proj;
  proj = new TLorentzVector(0,0,p_momentum, p_energy); //---------------PION BEAM
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

  TLorentzVector *prot1;
  prot1 = new TLorentzVector(0,0,0,0);
  TLorentzVector *prot2;
  prot2 = new TLorentzVector(0,0,0,0);
  TLorentzVector *invm_pepem;
  invm_pepem = new TLorentzVector(0,0,0,0);
  TLorentzVector *invm_pp;
  invm_pp = new TLorentzVector(0,0,0,0);
  
  TLorentzVector *mm;
  mm = new TLorentzVector(0,0,0,0);
  TLorentzVector *invm;
  invm = new TLorentzVector(0,0,0,0);
  TLorentzVector *invmCM;
  invmCM = new TLorentzVector(0,0,0,0);
  TLorentzVector *total;
  total = new TLorentzVector(0,0,0,0);
		   
  Double_t gloW; 

   TH1F *hepemInvM = new TH1F("hepemInvM","hepemInvM",111,0.0,0.555);
   hepemInvM->Sumw2();
   TH1F *hepemInvM_ps = new TH1F("hepemInvM_ps","hepemInvM_ps",111,0.0,0.555);
   hepemInvM_ps->Sumw2();
   TH1F *hepemInvM_div = new TH1F("hepemInvM_div","hepemInvM_div",111,0.0,0.555);
   hepemInvM_div->Sumw2();

   
   TH1F *hepemInvM1 = new TH1F("hepemInvM_nocut","hepemInvM_nocut",100,0.,0.6);
   TH1F *hepemMM = new TH1F("hepemMM","hepemMM",100,0.,1.2);
   TH1F *hOA = new TH1F("hOA","hOA",180, 0.0, 180);

   TH1F *hydil = new TH1F("hydil","hydil",100, -2.5, 2.5);
   TH1F *hptdil = new TH1F("hptdil","hptdil",100, 0, 0.5);
   TH1F *hydil_ps = new TH1F("hydil_ps","hydil_ps",100, -2.5, 2.5);
   TH1F *hptdil_ps = new TH1F("hptdil_ps","hptdil_ps",100, 0, 0.5);
   TH1F *hydil_div = new TH1F("hydil_div","hydil_div",100, -2.5, 2.5);
   TH1F *hptdil_div = new TH1F("hptdil_div","hptdil_div",100, 0, 0.5);

   TH1F *hEdil = new TH1F("hEdil","hEdil",100, 0, 1.);
   TH1F *hEdil1 = new TH1F("hEdil_nocut","hEdil_nocut",100, 0, 1.);

   TH2F *sigInvM = new TH2F("sigInvM","sigInvM",100,0.,0.4,100,0,3);
   TH2F *yInvM = new TH2F("yInvM","yInvM",100,0.,0.4,100,-2.5,2.5);
   TH2F *ptInvM = new TH2F("ptInvM","ptInvM",100,0.,0.4,100,0,0.5);
   
   TH2F *hM2M2 = new TH2F("hM2M2","hM2M2",100, 0., 2., 100, 0, 6.0);
   TH2F *hMpt = new TH2F("hMpt","hMpt",100, 0., 0.5, 100, 0, 0.6);
   TH2F *hMy = new TH2F("hMy","hMy",100, -2.5, 2.5, 100, 0, 0.6);
   TH3F *hMEdOm = new TH3F("hMEdOm","hMEdOm", 100 , 0, 0.6, 100,0.0, 0.6, 360, 0.0, 360.0);
   
   TH3F *epem_4pi_Wolf_wt = new TH3F("epem_4pi_Wolf_wt_Mypt","epem_4pi_Wolf_wt_Mypt", 111,0.0,0.555, 100,-2.5, 2.5, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_ps = new TH3F("epem_4pi_Wolf_ps_Mypt","epem_4pi_Wolf_ps_Mypt", 111,0.0,0.555, 100,-2.5, 2.5, 100, 0, 0.5);      

   TH2F *hM2M2_ps = new TH2F("hM2M2_ps","hM2M2_ps",100, 0., 2., 100, 0, 6.0);
   TH2F *hMpt_ps = new TH2F("hMpt_ps","hMpt_ps",100, 0., 0.5, 100, 0, 0.6);
   TH2F *hMy_ps = new TH2F("hMy_ps","hMy_ps",100, -2.5, 2.5, 100, 0, 0.6);
   TH3F *hMEdOm_ps = new TH3F("hMEdOm_ps","hMEdOm_ps", 100 , 0, 0.6, 100,0.0, 0.6, 360, 0.0, 360.0);

   TH2F *hM2M2_div = new TH2F("hM2M2_div","hM2M2_div",100, 0., 2., 100, 0, 6.0);
   TH2F *hMpt_div = new TH2F("hMpt_div","hMpt_div",100, 0., 0.5, 100, 0, 0.6);
   TH2F *hMy_div = new TH2F("hMy_div","hMy_div",100, -2.5, 2.5, 100, 0, 0.6);
   TH3F *hMEdOm_div = new TH3F("hMEdOm_div","hMEdOm_div", 100 , 0, 0.6, 100,0.0, 0.6, 360, 0.0, 360.0);

   
   float sig[8]={0.};
   float nb[8]={0.};
   TGraph *gr1;

  float nbEv=1000000.;

  
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
	 
  TVector3 m3,vert, p1, p2, p3, p4;
  TLorentzVector *lv = new TLorentzVector;
  TLorentzVector *Theta_CM;
  TLorentzVector *lv1 = new TLorentzVector;
  TLorentzVector *lv2 = new TLorentzVector;	
  TLorentzVector *lv3 = new TLorentzVector;
  TLorentzVector *lv4 = new TLorentzVector;
  
  //m3.SetXYZ(Particles_fP_fX[i],Particles_fP_fY[i],Particles_fP_fZ[i]);
  //vert.SetXYZ(Particles_fV_fX[i],Particles_fV_fY[i],Particles_fV_fZ[i]);

  lv1->SetPxPyPzE(a2,a3,a4,a1);	 
  lv2->SetPxPyPzE(b2,b3,b4,b1);	 
  lv3->SetPxPyPzE(c2,c3,c4,c1);	 
  lv4->SetPxPyPzE(d2,d3,d4,d1);	 

    

    p1.SetXYZ(lv1->Px(),lv1->Py(),lv1->Pz());
    prot1->SetVectM(p1, mpi);

    p2.SetXYZ(lv2->Px(),lv2->Py(),lv2->Pz());
    prot2->SetVectM(p2, mp);
  
    p3.SetXYZ(lv3->Px(),lv3->Py(),lv3->Pz());
    ep->SetVectM(p3, me);

    p4.SetXYZ(lv4->Px(),lv4->Py(),lv4->Pz());
    em->SetVectM(p4, me);

    *invm_pp=*prot1 + *prot2;
    *invm_pepem=*prot1 + *ep + *em;
    
    *mm=*beam-*ep-*em;
    *invm = *ep + *em;
    
    float Invmass = invm->M();
    float InvEn = invm->E();

    float InvM_pp = invm_pp->M2();
    float InvM_pepem = invm_pepem->M2();

    //cout<<InvM_pp<<" "<<InvM_pepem<<endl;
    //float Ep=ep->E();
    //float Em=em->E();
    float en1=InvEn*sqrt(InvEn*InvEn-Invmass*Invmass);
    //cout<<InvEn<<" "<<Ep<<" "<<Em<<" "<<Ep+Em<<endl;
      
    float InvMom = invm->P();
    float pt = invm->Pt();
    float y = invm->Rapidity();
    float y_ep = ep->Rapidity();
    float y_em = em->Rapidity();
    float Invtheta = invm->Theta()*TMath::RadToDeg();
    float Invphi = invm->Phi()*TMath::RadToDeg();
    float dOmega = TMath::Cos(Invtheta)*Invtheta*Invphi;
    //    cout<<Invtheta<<endl;
    float pt_ep = ep->Pt();
    float pt_em = em->Pt();
    float missmass = mm->M();


    ep->Boost(-(*beam).BoostVector() );
    em->Boost(-(*beam).BoostVector() );
    *invmCM = *ep + *em;
    
    float InvmassCM = invmCM->M();
    float InvEnCM = invmCM->E();
    float yCM = invmCM->Rapidity();
    float ptCM = invmCM->Pt();
    float pCM = invmCM->P();
    
    //cout<<InvEn<<" "<<InvEnCM<<endl;
    //cout<<"y: "<<y<<" "<<yCM<<endl;
    //cout<<"pt: "<<pt<<" "<<ptCM<<endl;


    
    Double_t Opang = openingangle(*ep,*em);
    Double_t Opang_deg = Opang*TMath::RadToDeg();

    float Energy = sqrt(pow(Invmass,2)+pow(pt,2))*TMath::CosH(y);

    float s=sqrts*sqrts;
    float s2=s+(Invmass*Invmass)-(2.*InvEnCM*sqrts);
    //float s2=s+(Invmass*Invmass)-(2.*Energy*sqrts);

    //float s2=s+(Invmass*Invmass)-(2.*en1*sqrts);

    float R2s=sqrt(1.-((mp+mpi)*(mp+mpi)/s));
    float R2s2=sqrt(1.- ( (mp+mpi)*(mp+mpi)/s2));
    //float R=sqrt((s*(s2-4*mp*mp))/(s2*(s-4*mp*mp)));
    float sigma_quasiel= ((s-(mp+mpi)*(mp+mpi))/(2.*mpi*mpi))*sigma_el;
    //float dsigma=A*(sigma_quasiel/(Invmass*InvEn*InvEn))*(R2s2/R2s)*2*pt;
    //float dsigma=A*(sigma_quasiel/(Invmass*en1*en1))*(R2s2/R2s);
    //----------float dsigma=A*(sigma_quasiel*InvMom/(Invmass*InvEn*InvEn))*(R2s2/R2s);
    float dsigma=A*pCM*(sigma_quasiel/(InvmassCM*InvEnCM*InvEnCM))*(R2s2/R2s);

    //cout<<"::: "<<R2s2<<"  "<<((mpi+mp)*(mpi+mp)/s2)<<endl;

    
    if (((mp+mpi)*(mp+mpi)/s2)>1)continue;
    // cout<<(R2s2/R2s)<<endl;
    // if (TMath::IsNaN(dsigma))continue;

    hepemInvM_ps->Fill(Invmass,(1./nbEv));
    hydil_ps->Fill(y,1./nbEv);
    hptdil_ps->Fill(pt,1./nbEv);
 
    hepemInvM->Fill(Invmass,dsigma*(1./nbEv));
    hydil->Fill(y,dsigma*(1./nbEv));
    hptdil->Fill(pt,dsigma*(1./nbEv));
    //hEdil->Fill(InvEn);
    

    hMpt->Fill(pt,Invmass,dsigma*(1./nbEv));
    hMy->Fill(y,Invmass,dsigma*(1./nbEv));
    hM2M2->Fill(InvM_pepem,InvM_pp,dsigma*(1./nbEv));
    hMEdOm->Fill(Invmass,InvEnCM,dOmega,dsigma*(1./nbEv));       

    hMpt_ps->Fill(pt,Invmass,(1./nbEv));
    hMy_ps->Fill(y,Invmass,(1./nbEv));
    hM2M2_ps->Fill(InvM_pepem,InvM_pp,(1./nbEv));
    hMEdOm_ps->Fill(Invmass,InvEnCM,dOmega,dsigma*(1./nbEv));
    
    epem_4pi_Wolf_wt->Fill(Invmass,InvEnCM,dOmega,dsigma*(1./nbEv));       
    epem_4pi_Wolf_ps->Fill(Invmass,InvEnCM,dOmega,(1./nbEv));       

 }


 hepemInvM_div->Divide(hepemInvM,hepemInvM_ps,1.,1.);
 hydil_div->Divide(hydil,hydil_ps,1.,1.);
 hptdil_div->Divide(hptdil,hptdil_ps,1.,1.);
    

 hMpt_div->Divide(hMpt,hMpt_ps,1.,1.);
 hMy_div->Divide(hMy,hMy_ps,1.,1.);
 hM2M2_div->Divide(hM2M2,hM2M2_ps,1.,1.);
 hMEdOm_div->Divide(hMEdOm,hMEdOm_ps,1.,1.);


 //  TGraph *gr=new TGraph(8,invM,sigfin);
 // gr->SetMarkerStyle(21);
 // gr->SetMarkerSize(1.);
 hMEdOm->Draw("ISO");
 TH1* hMEdOm_px = hMEdOm->ProjectionX("hMEdOm_px");
 hMEdOm_px->Draw();
 TH1* hMEdOm_ps_px = hMEdOm_ps->ProjectionX("hMEdOm_ps_px");
 hMEdOm_ps_px->Draw();

 //TCanvas *can1=new TCanvas("brems_theory1","brems_theory1");
 //can1->cd();
 //gPad->SetLogy();
 //gr->Draw("al");

 TH1* epem_4pi_Wolf_wt_px = epem_4pi_Wolf_wt->ProjectionX("epem_4pi_Wolf_wt_px");
 epem_4pi_Wolf_wt_px->Draw();
 TH1* epem_4pi_Wolf_ps_px = epem_4pi_Wolf_ps->ProjectionX("epem_4pi_Wolf_ps_px");
 epem_4pi_Wolf_ps_px->Draw();

 
 float xWolf[]={0.05565393,0.06698674,0.07376166,0.08132172,0.089637786,0.09870986,0.10778921,0.11763914,0.12673306, 0.1350564, 0.14492819, 0.15325882, 0.16236003, 0.17146122, 0.18131845, 0.18966363,0.19799428,0.20861477, 0.21847199, 0.2275659,0.23819369, 0.2480509, 0.25944197, 0.2692846,0.28066114,0.29126707,0.30188027,0.31171563,0.32003897,0.32987434, 0.3389537, 0.3502938, 0.3578247,0.3653848, 0.37217426,0.3804612,0.38722154,0.3939819,0.39848882};

 float yWolf[]={0.002029871, 0.0015116982, 0.0011904291, 9.900433E-4,8.0835074E-4, 6.4794865E-4, 5.289451E-4, 4.396775E-4,3.7227408E-4, 3.0955495E-4, 2.7179986E-4, 2.301725E-4, 1.9847752E-4, 1.7114694E-4, 1.4488454E-4,  1.2725792E-4,1.0776782E-4, 9.121501E-5, 7.721812E-5, 6.538044E-5, 5.635784E-5, 4.770976E-5, 4.111858E-5, 3.3560813E-5, 2.7887192E-5, 2.2757444E-5,1.891348E-5, 1.5157812E-5, 1.26040895E-5, 1.0101283E-5, 8.246061E-6, 6.2542103E-6,  4.835106E-6, 4.0212085E-6, 3.2843816E-6, 2.49277E-6,  1.8926141E-6, 1.436951E-6, 1.195901E-6};


 
 TGraph *gr2 = new TGraph(36,xWolf,yWolf);
 gr2->SetMarkerStyle(21);
 gr2->SetMarkerSize(1.);
 gr2->SetMarkerColor(2);
  
  
 //TCanvas *can2=new TCanvas("brems_th_wolf_en1","brems_th_wolf_en1");
 TCanvas *can2=new TCanvas("brems_th_wolf_En","brems_th_wolf_En");

 can2->cd();
 gPad->SetLogy();
 //hepemInvM->Scale(1./(1000000.*6));
 //hepemInvM->Scale(1./(1000000.));

 //hepemInvM->Scale(1./(1000000.));
 //hepemInvM->Scale(1./hepemInvM->GetEntries());   
 normalize(hepemInvM);
 normalize(hepemInvM_ps);
 hepemInvM->Draw();
 gr2->Draw("psame");


 TCanvas *can3=new TCanvas("brems_2D","brems_2D",900,700);
 can3->Divide(3,4);

 can3->cd(1);
 hMpt_ps->Draw("colz");
 can3->cd(2);
 hMpt->Draw("colz");
 can3->cd(3);
 hMpt_div->Draw("colz");

 can3->cd(4);
 hMy_ps->Draw("colz");
 can3->cd(5);
 hMy->Draw("colz");
 can3->cd(6);
 hMy_div->Draw("colz");

 can3->cd(7);
 hM2M2_ps->Draw("colz");
 can3->cd(8);
 hM2M2->Draw("colz");
 can3->cd(9);
 hM2M2_div->Draw("colz");

 can3->cd(10);
 hMEdOm_ps->SetFillColor(kCyan);
 hMEdOm_ps->Draw("ISO");
 can3->cd(11);
 hMEdOm->SetFillColor(kCyan);
 hMEdOm->Draw("ISO");
 can3->cd(12);
 hMEdOm_div->SetFillColor(kCyan);
 hMEdOm_div->Draw("ISO");



 TFile *output = new TFile("Wolf_caln_PimP_Energy_GiBUU_CHECK.root","RECREATE");
 output->cd();

 hMEdOm_px->Write();
 hMEdOm_ps_px->Write(); 
 TH1F *Wolf_Brem = (TH1F*)hMEdOm_px->Clone();
 Wolf_Brem->Divide(hMEdOm_ps_px);
 Wolf_Brem->SetName("Wolf_Brem");
 Wolf_Brem->Draw("same,hist");
 Wolf_Brem->Write();


 epem_4pi_Wolf_wt_px->Write();
 epem_4pi_Wolf_ps_px->Write();
 TH1F *Wolf_Brem_2 = (TH1F*)epem_4pi_Wolf_wt_px->Clone();
 Wolf_Brem_2->Divide(epem_4pi_Wolf_ps_px);
 Wolf_Brem_2->SetName("Wolf_Brem2");
 Wolf_Brem_2->Draw("same,hist");
 Wolf_Brem_2->Write();


 hepemInvM->Write();
 hepemInvM_ps->Write(); 
 TH1F *Wolf_Brem_3 = (TH1F*)hepemInvM->Clone();
 Wolf_Brem_3->Divide(hepemInvM_ps);
 Wolf_Brem_3->SetName("Wolf_Brem_GiBUU");
 Wolf_Brem_3->Draw("same,hist");
 Wolf_Brem_3->Write();

 output->Write();
 output->Close();
 
 
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



 
