#include "TPaveText.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1.h"
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

void wolf_pimp_brem_dSd3p(){
  ifstream f1a ;
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PPim_Bremsstrahlung_P690MeV.evt");
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/690/PHASESPACE/pion690_brems_4.evt");
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/690/PHASESPACE/pion690_BREM_22M.evt");
  f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/690/ANGULAR/pion690_BREM_ANG_22M.evt");
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/690/PimN_PimN_Brem_PhS/pimn_690_pimn_Brem_22M.evt");
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/690/PimP_Pi0N_Brem_PhS/pimp_690_pi0n_Brem_22M.evt");
  
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_PPim_brems_new/p_ee_0.evt");

  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/CARBON/pimC_Brem_PHASESPACE/Pim12C_Brem_PHASESPACE_22_FULL.evt");
 const float PI= TMath::Pi();

 double p_momentum =  0.690;//GeV/c
 const float mpi=0.13957;
 const float mp=0.93827;
 //const float mpi=0.93827;
 //const float mp=0.93956;
 const float me=0.510998*0.001;

 double p_energy = sqrt(p_momentum*p_momentum + mpi*mpi);

 const float alpha = 1./137.;
 const float sqrts=1.5;
 const float sigma_el=18.130;        //mb    pi- p --> pi- p e+e-
 //const float sigma_el=9.77;       //mb    pi- p --> pi0 n  e+e-
 //const float sigma_el=12.1;       //mb    pi- n --> pi- n  e+e- 
 //const float sigma_el=43.10;      //mb 
 const float A=alpha*alpha/(6.*PI*PI);
 const float A_SCH=2.*alpha*alpha/(3.*PI*PI);
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

   TH1F *hepemInvM = new TH1F("hepemInvM","hepemInvM",100,0.,0.6);
   hepemInvM->Sumw2();
   TH1F *hepemInvMCM = new TH1F("hepemInvMCM","hepemInvMCM",100,0.,0.6);
   hepemInvMCM->Sumw2();

   TH1F *hepemInvM_ps = new TH1F("hepemInvM_ps","hepemInvM_ps",100,0.,0.6);
   hepemInvM_ps->Sumw2();
   TH1F *hepemInvM_div = new TH1F("hepemInvM_div","hepemInvM_div",100,0.,0.6);
   hepemInvM_div->Sumw2();

   
   TH1F *hepemInvM1 = new TH1F("hepemInvM_nocut","hepemInvM_nocut",100,0.,0.6);
   TH1F *hepemMM = new TH1F("hepemMM","hepemMM",100,0.,1.2);
   TH1F *hOA = new TH1F("hOA","hOA",180, 0.0, 180);

   TH1F *hydil = new TH1F("hydil","hydil",200, -5.0, 5.0);
   TH1F *hptdil = new TH1F("hptdil","hptdil",100, 0, 0.5);
   TH1F *hydilCM = new TH1F("hydilCM","hydilCM",200, -5.0, 5.0);
   TH1F *hptdilCM = new TH1F("hptdilCM","hptdilCM",100, 0, 0.5);
   TH1F *hydil_ps = new TH1F("hydil_ps","hydil_ps",100, -3.0,4.0);
   TH1F *hptdil_ps = new TH1F("hptdil_ps","hptdil_ps",100, 0, 0.5);
   TH1F *hydil_div = new TH1F("hydil_div","hydil_div",100, -3.0,4.0);
   TH1F *hptdil_div = new TH1F("hptdil_div","hptdil_div",100, 0, 0.5);

   TH1F *hEdil = new TH1F("hEdil","hEdil",100, 0, 1.);
   TH1F *hEdil1 = new TH1F("hEdil_nocut","hEdil_nocut",100, 0, 1.);

   TH2F *sigInvM = new TH2F("sigInvM","sigInvM",100,0.,0.4,100,0,3);
   TH2F *yInvM = new TH2F("yInvM","yInvM",100,0.,0.4,100,-3.0,4.0);
   TH2F *ptInvM = new TH2F("ptInvM","ptInvM",100,0.,0.4,100,0,0.5);
   
   TH2F *hM2M2 = new TH2F("hM2M2","hM2M2",100, 0., 2., 100, 0, 6.0);
   TH2F *hMpt = new TH2F("hMpt","hMpt",100, 0., 0.5, 100, 0, 0.6);
   TH2F *hMy = new TH2F("hMy","hMy",100,-3.0,4.0, 100, 0, 0.6);
   TH3F *hMypt = new TH3F("hMypt","hMypt", 100 , 0, 0.6, 100,-3.0,4.0, 100, 0, 0.5);
   
   TH3F *epem_4pi_Wolf_wt = new TH3F("epem_4pi_Wolf_wt_Mypt","epem_4pi_Wolf_wt_Mypt", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_M_wt = new TH3F("epem_4pi_Wolf_M_wt_Mypt","epem_4pi_Wolf_M_wt_Mypt", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_Ecm_wt = new TH3F("epem_4pi_Wolf_Ecm_wt_Mypt","epem_4pi_Wolf_Ecm_wt_Mypt", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_ps = new TH3F("epem_4pi_Wolf_ps_Mypt","epem_4pi_Wolf_ps_Mypt", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_wt_ONLY = new TH3F("epem_4pi_Wolf_wt_Mypt_ONLY","epem_4pi_Wolf_wt_Mypt_ONLY", 200,0.0,1.0, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_wt_DIV = new TH3F("epem_4pi_Wolf_wt_Mypt_DIV","epem_4pi_Wolf_wt_Mypt_DIV", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_M_wt_DIV = new TH3F("epem_4pi_Wolf_M_wt_Mypt_DIV","epem_4pi_Wolf_M_wt_Mypt_DIV", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_Ecm_wt_DIV = new TH3F("epem_4pi_Wolf_Ecm_wt_Mypt_DIV","epem_4pi_Wolf_Ecm_wt_Mypt_DIV", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);      
   
   TH3F *epem_4pi_Wolf_wt1 = new TH3F("epem_4pi_Wolf_wt_Mypt1","epem_4pi_Wolf_wt_Mypt1", 600,0.0,0.600, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_ps1 = new TH3F("epem_4pi_Wolf_ps_Mypt1","epem_4pi_Wolf_ps_Mypt1", 600,0.0,0.600, 100,-3.0,4.0, 100, 0, 0.5);      
/*
   TH3F *epem_4pi_Wolf_wt_special = new TH3F("epem_4pi_Wolf_wt_Mypt_special","epem_4pi_Wolf_wt_Mypt_special", 200,0.0,1.0, 100,-5.0,5.0, 200, 0, 1.0);      
   TH3F *epem_4pi_Wolf_ps_special = new TH3F("epem_4pi_Wolf_ps_Mypt_special","epem_4pi_Wolf_ps_Mypt_special", 200,0.0,1.0, 100,-5.0,5.0, 200, 0, 1.0);      
   TH3F *epem_4pi_Wolf_wt_special_DIV = new TH3F("epem_4pi_Wolf_wt_Mypt_special_DIV","epem_4pi_Wolf_wt_Mypt_special_DIV", 111,0.0,0.555, 100,-3.0,4.0, 100, 0, 0.5);
*/
   TH3F *epem_4pi_Wolf_wt_special = new TH3F("epem_4pi_Wolf_wt_Mypt_special","epem_4pi_Wolf_wt_Mypt_special", 100 , 0, 0.6, 100,-3.0,4.0, 100, 0, 0.5);
   TH3F *epem_4pi_Wolf_ps_special = new TH3F("epem_4pi_Wolf_ps_Mypt_special","epem_4pi_Wolf_ps_Mypt_special", 100 , 0, 0.6, 100,-3.0,4.0, 100, 0, 0.5);      
   TH3F *epem_4pi_Wolf_wt_special_DIV = new TH3F("epem_4pi_Wolf_wt_Mypt_special_DIV","epem_4pi_Wolf_wt_Mypt_special_DIV", 100 , 0, 0.6, 100,-3.0,4.0, 100, 0, 0.5);
   
   
   TH2F *hM2M2_ps = new TH2F("hM2M2_ps","hM2M2_ps",100, 0., 2., 100, 0, 6.0);
   TH2F *hMpt_ps = new TH2F("hMpt_ps","hMpt_ps",100, 0., 0.5, 100, 0, 0.6);
   TH2F *hMy_ps = new TH2F("hMy_ps","hMy_ps",100,-3.0,4.0, 100, 0, 0.6);
   TH3F *hMypt_ps = new TH3F("hMypt_ps","hMypt_ps", 100 , 0, 0.6, 100,-3.0,4.0, 100, 0, 0.5);      

   TH2F *hM2M2_div = new TH2F("hM2M2_div","hM2M2_div",100, 0., 2., 100, 0, 6.0);
   TH2F *hMpt_div = new TH2F("hMpt_div","hMpt_div",100, 0., 0.5, 100, 0, 0.6);
   TH2F *hMy_div = new TH2F("hMy_div","hMy_div",100, -3.0,4.0, 100, 0, 0.6);
   TH3F *hMypt_div = new TH3F("hMypt_div","hMypt_div", 100 , 0, 0.6, 100,-3.0,4.0, 100, 0, 0.5);      

   TH1F *PhS1 = new TH1F("PS_OG","PS_OG",150,0.,1.5);PhS1->Sumw2();
   TH1F *PhS1_M= new TH1F("PS_InvM","PS_InvM",150,0.,1.5);PhS1_M->Sumw2();
   TH1F *PhS1_Ecm= new TH1F("PS_Ecm","PS_Ecm",150,0.,1.5);PhS1_Ecm->Sumw2();
   TH1F *DS_inv = new TH1F("DS_InvM","DS_Inv",100,0.,1.0); DS_inv->Sumw2();
   TH1F *DS_M = new TH1F("DS_M","DS_M",100,0.,1.0); DS_M->Sumw2();
   TH1F *DS_Ecm = new TH1F("DS_Ecm","DS_Ecm",100,0.,1.0); DS_Ecm->Sumw2();

   TH1F *hepemInvM_PLUTO1 = new TH1F("hepemInvM_PLUTO1","hepemInvM_PLUTO1",111,0.0,0.555);
   TH1F *DiffSec_saf = new TH1F("DiffSec_saf","DiffSec_saf",200,0.,1.0); DiffSec_saf->Sumw2();
   
   float sig[8]={0.};
   float nb[8]={0.};
   TGraph *gr1;

   float nbEv=2200000.;
   float nbEv1=1645986.;

  //ofstream f1;
  //f1.open ("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/PPIM_PHASESPACE/pion690_brems_wts_4.evt");
  //f1.open ("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/PPIM_PHASESPACE/pion690_BREM_22M_WEIGHTs_WOLFs.evt");
  
   float ii = 0;  
while (!f1a.eof()){


  float p1a,wa,w1a,w2a,w3a;
  float a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13;
  float b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13;
  float c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13;
  float d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13;

  f1a>>p1a>>wa>>w1a>>w2a>>w3a;
  //f1a>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12>>a13;
  f1a>>a1>>a2>>a3>>a4>>a5>>a6>>a7>>a8>>a9>>a10>>a11>>a12>>a13;
  f1a>>b1>>b2>>b3>>b4>>b5>>b6>>b7>>b8>>b9>>b10>>b11>>b12>>b13;
  f1a>>c1>>c2>>c3>>c4>>c5>>c6>>c7>>c8>>c9>>c10>>c11>>c12>>c13;
  f1a>>d1>>d2>>d3>>d4>>d5>>d6>>d7>>d8>>d9>>d10>>d11>>d12>>d13;

	 
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
    float TH = invm->Theta() * TMath::RadToDeg();
    float theta = TMath::Sin(TH);
    float phi = invm->Phi();
    float dOmega = 2*PI*theta*TH*phi;
    //cout<<InvEn<<" "<<InvEnCM<<endl;
    cout<<"theta: "<<theta<<"   "<<"TH: "<<TH<<"   "<<"phi :"<<phi<<"   "<<"Omega : "<<dOmega<<endl;
    //cout<<"pt: "<<pt<<" "<<ptCM<<endl;


    
    Double_t Opang = openingangle(*ep,*em);
    Double_t Opang_deg = Opang*TMath::RadToDeg();

    float Energy = sqrt(pow(Invmass,2)+pow(pt,2))*TMath::CosH(y);

    float s=sqrts*sqrts;
    
    float s2_Ecm=s+(Invmass*Invmass)-(2.*InvEnCM*sqrts);
    float s2_M=s+(Invmass*Invmass)-(2.*Invmass*sqrts);
    float s2=s+(Invmass*Invmass)-(2.*Energy*sqrts);

    //float s2=s+(Invmass*Invmass)-(2.*en1*sqrts);
    float Qcm = sqrt((pow((s+(mpi*mpi)-(mp*mp)),2)/(4*s))-(mpi*mpi));
    float Qcm_s2 = sqrt((pow((s2+(mpi*mpi)-(mp*mp)),2)/(4*s2))-(mpi*mpi));

    float Qcm1 = sqrt(pow((s+(mpi*mpi)-(mp*mp)),2)/((4*s)-(mpi*mpi)));
    float Qcm_s21 = sqrt(pow((s2+(mpi*mpi)-(mp*mp)),2)/((4*s2)-(mpi*mpi)));

    
    //float Qcm = sqrt((s-pow((mp+mpi),2))*(s-pow((mp-mpi),2)))/(2*sqrts);
    //float Qcm_s2 = sqrt((s2-pow((mp+mpi),2))*(s2-pow((mp-mpi),2)))/(2*sqrts);

    //float Qcm_s2_M = sqrt((s2_M-pow((mp+mpi),2))*(s2_M-pow((mp-mpi),2)))/(2*sqrt(s2_M));
    //float Qcm_s2_Ecm = sqrt((s2_Ecm-pow((mp+mpi),2))*(s2_Ecm-pow((mp-mpi),2)))/(2*sqrt(s2_Ecm));
    
    float Qcm_s2_M = sqrt((pow((s2_M+(mpi*mpi)-(mp*mp)),2)/(4*s2_M))-(mpi*mpi));
    float Qcm_s2_Ecm = sqrt((pow((s2_Ecm+(mpi*mpi)-(mp*mp)),2)/(4*s2_Ecm))-(mpi*mpi));
    
    //float R2s=sqrt(1.-(((mp+mpi)*(mp+mpi))/s));
    float R2s1=(2*Qcm)/sqrts;
    float R2s21= (2*Qcm_s2)/sqrt(s2);
    float R2s2_M= (2*Qcm_s2_M)/sqrt(s2_M);
    float R2s2_Ecm= (2*Qcm_s2_Ecm)/sqrt(s2_Ecm);

    float R2s=sqrt(1.-(((mp+mpi)*(mp+mpi))/s));
    float R2s2=sqrt(1.- (((mp+mpi)*(mp+mpi))/s2));
    //float R2s2_M=sqrt(1.-(((mp+mpi)*(mp+mpi))/s2_M));
    //float R2s2_Ecm=sqrt(1.-(((mp+mpi)*(mp+mpi))/s2_Ecm));
    float PS1=(R2s2/R2s);          
    float PS1_M=(R2s2_M/R2s);          
    float PS1_Ecm=(R2s2_Ecm/R2s);          

    //float R=sqrt((s*(s2-4*mp*mp))/(s2*(s-4*mp*mp)));
    float sigma_quasiel= ((s-((mp+mpi)*(mp+mpi)))/(2.*mpi*mpi))*sigma_el;

    float sigma_quasiel1= ((Qcm*Qcm)/(mpi*mpi))*2*sigma_el;
    //float sigma_quasiel= ((s-((mpi+mp)*(mpi+mp)))/((mpi*mpi)+(mp*mp)))*18.13;
    //cout<<"W0lf's sig :  "<<sigma_quasiel<<"GiBUU sig :  "<<sigma_quasiel1<<endl;    
    //float dsigma=A*(sigma_quasiel/(Invmass*InvEnCM*InvEnCM))*(R2s2/R2s)*2*pt;
    //float dsigma=A*(sigma_quasiel/(Invmass*en1*en1))*(R2s2/R2s);
    float dsigma=A*(sigma_quasiel1/(Invmass*Energy*Energy))*2*pt;
    //float dsigma=A*(sigma_quasiel1/(Invmass*Energy*Energy))*(R2s21/R2s1)*2*pt;
    float dsigma_M=A*(sigma_quasiel1/(Invmass*Energy*Energy))*(R2s2_M/R2s1)*2*pt;
    float dsigma_Ecm=A*(sigma_quasiel1/(Invmass*Energy*Energy))*(R2s2_Ecm/R2s1)*2*pt;

    //float weights = (1/dsigma);

/////////////////////////////////---------------SCHAFER--------------------////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //if (((mp+mpi)*(mp+mpi)/s2)>1)continue;
    // cout<<(R2s2/R2s)<<endl;
    if (TMath::IsNaN(dsigma))continue;
    std::cout.precision(20);
    //cout<<p1a<<endl;
    ii=ii+1;
    cout << ii << endl;
    //cout<<p1a<<"::: "<<(R2s2/R2s)<<":::::"<<ii<<endl;

    hepemInvM_ps->Fill(Invmass);
    hydil_ps->Fill(y);
    hptdil_ps->Fill(pt);
 
    hepemInvM->Fill(Invmass,dsigma);
    hydil->Fill(y,dsigma);
    hptdil->Fill(pt,dsigma);

    hepemInvMCM->Fill(InvmassCM,dsigma);
    hydilCM->Fill(yCM,dsigma);
    hptdilCM->Fill(ptCM,dsigma);

    PhS1->Fill(PS1);
    PhS1_M->Fill(PS1_M);
    PhS1_Ecm->Fill(PS1_Ecm);
    //hEdil->Fill(InvEn);
    DS_inv->Fill(dsigma);
    DS_M->Fill(dsigma_M);
    DS_Ecm->Fill(dsigma_Ecm);

    hMpt->Fill(pt,Invmass,dsigma);
    hMy->Fill(y,Invmass,dsigma);
    hM2M2->Fill(InvM_pepem,InvM_pp,dsigma);
    hMypt->Fill(Invmass,y,pt,dsigma);       

    hMpt_ps->Fill(pt,Invmass);
    hMy_ps->Fill(y,Invmass);
    hM2M2_ps->Fill(InvM_pepem,InvM_pp);
    hMypt_ps->Fill(Invmass,y,pt);       

    epem_4pi_Wolf_wt->Fill(Invmass,y,pt,dsigma);       
    epem_4pi_Wolf_M_wt->Fill(Invmass,y,pt,dsigma_M);
    epem_4pi_Wolf_Ecm_wt->Fill(Invmass,y,pt,dsigma_Ecm);
    epem_4pi_Wolf_ps->Fill(Invmass,y,pt);       
    epem_4pi_Wolf_wt_ONLY->Fill(dsigma,y,pt);       

    epem_4pi_Wolf_wt1->Fill(Invmass,y,pt,dsigma);       
    epem_4pi_Wolf_ps1->Fill(Invmass,y,pt);       

    epem_4pi_Wolf_wt_special->Fill(Invmass,y,pt,dsigma);       
    epem_4pi_Wolf_ps_special->Fill(Invmass,y,pt);       
 //////////////////////////////////////////////////////////////////////////////////////////////////////   
 //////////////////////////////////////////////////////////////////////////////////////////////////////   
/*
    cout<<p1a<<" "<<wa<<" "<<w1a<<" "<<w2a<<" "<<w3a<<endl;
    cout<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a13*dsigma<<endl;
    cout<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<" "<<b11<<" "<<b12<<" "<<b13*dsigma<<endl;
    cout<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<" "<<c11<<" "<<c12<<" "<<c13*dsigma<<endl;
    cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<" "<<d7<<" "<<d8<<" "<<d9<<" "<<d10<<" "<<d11<<" "<<d12<<" "<<d13*dsigma<<endl;

    f1<<p1a<<" "<<wa<<" "<<w1a<<" "<<w2a<<" "<<w3a<<endl;
    f1<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<dsigma<<"\n";
    f1<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<" "<<b11<<" "<<b12<<" "<<dsigma<<"\n";
    f1<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<" "<<c11<<" "<<c12<<" "<<dsigma<<"\n";
    f1<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<" "<<d7<<" "<<d8<<" "<<d9<<" "<<d10<<" "<<d11<<" "<<d12<<" "<<dsigma<<"\n";
*/
 //////////////////////////////////////////////////////////////////////////////////////////////////////   
 //////////////////////////////////////////////////////////////////////////////////////////////////////   
 }
//DS_invM->Draw();
 TCanvas *can3D=new TCanvas("brems_can3D","brems_can3D");
 can3D->cd();
 hepemInvM_div->Divide(hepemInvM,hepemInvM_ps,1.,1.);
 hydil_div->Divide(hydil,hydil_ps,1.,1.);
 hptdil_div->Divide(hptdil,hptdil_ps,1.,1.);

//------------------------------------------------------------------------------
 epem_4pi_Wolf_wt_DIV->Divide(epem_4pi_Wolf_wt,epem_4pi_Wolf_ps,1.,1.);       
 epem_4pi_Wolf_wt_DIV->Draw("ISO");       
 TH1* epem_4pi_Wolf_wt_DIV_px = epem_4pi_Wolf_wt_DIV->ProjectionX("epem_4pi_Wolf_wt_DIV_px");
 epem_4pi_Wolf_wt_DIV_px->Draw();
//------------------------------------------------------------------------------
 epem_4pi_Wolf_M_wt_DIV->Divide(epem_4pi_Wolf_M_wt,epem_4pi_Wolf_ps,1.,1.);       
 epem_4pi_Wolf_M_wt_DIV->Draw("ISO");       
 TH1* epem_4pi_Wolf_M_wt_DIV_px = epem_4pi_Wolf_M_wt_DIV->ProjectionX("epem_4pi_Wolf_M_wt_DIV_px");
 epem_4pi_Wolf_M_wt_DIV_px->Draw();
//------------------------------------------------------------------------------
 epem_4pi_Wolf_Ecm_wt_DIV->Divide(epem_4pi_Wolf_Ecm_wt,epem_4pi_Wolf_ps,1.,1.);       
 epem_4pi_Wolf_Ecm_wt_DIV->Draw("ISO");       
 TH1* epem_4pi_Wolf_Ecm_wt_DIV_px = epem_4pi_Wolf_Ecm_wt_DIV->ProjectionX("epem_4pi_Wolf_Ecm_wt_DIV_px");
 epem_4pi_Wolf_Ecm_wt_DIV_px->Draw();
//------------------------------------------------------------------------------

 //epem_4pi_Wolf_wt_DIV = (TH3F*)epem_4pi_Wolf_wt->Clone();   
 //epem_4pi_Wolf_wt_DIV->SetName("epem_4pi_Wolf_wt_DIV");
 //epem_4pi_Wolf_wt_DIV->Draw("ISO");

 //epem_4pi_Wolf_wt_special_DIV = (TH3F*)epem_4pi_Wolf_wt_special->Clone();   
 //epem_4pi_Wolf_wt_special_DIV->Divide(epem_4pi_Wolf_ps_special); 
 //epem_4pi_Wolf_wt_special_DIV->SetName("epem_4pi_Wolf_wt_special_DIV");
 epem_4pi_Wolf_wt_special_DIV->Divide(epem_4pi_Wolf_wt_special,epem_4pi_Wolf_ps_special,1.,1.);       
 epem_4pi_Wolf_wt_special_DIV->Draw("ISO");

 
 hMpt_div->Divide(hMpt,hMpt_ps,1.,1.);
 hMy_div->Divide(hMy,hMy_ps,1.,1.);
 hM2M2_div->Divide(hM2M2,hM2M2_ps,1.,1.);
 hMypt_div->Divide(hMypt,hMypt_ps,1.,1.);
 TH1* hMypt_div_px = hMypt_div->ProjectionX("hMypt_div_px");
 hMypt_div_px->Draw();
 
 TCanvas *can3D1=new TCanvas("brems_can3D1","brems_can3D1");
 can3D1->cd();
 //  TGraph *gr=new TGraph(8,invM,sigfin);
 // gr->SetMarkerStyle(21);
 // gr->SetMarkerSize(1.);
 hMypt->Draw("ISO");
 
 TH1* hMypt_px = hMypt->ProjectionX("hMypt_px");
 hMypt_px->Draw();
 TH1* hMypt_ps_px = hMypt_ps->ProjectionX("hMypt_ps_px");
 hMypt_ps_px->Draw();

 //TCanvas *can1=new TCanvas("brems_theory1","brems_theory1");
 //can1->cd();
 //gPad->SetLogy();
 //gr->Draw("al");

 TH1* epem_4pi_Wolf_wt_px = epem_4pi_Wolf_wt->ProjectionX("epem_4pi_Wolf_wt_px");
 epem_4pi_Wolf_wt_px->Draw();
 TH1* epem_4pi_Wolf_ps_px = epem_4pi_Wolf_ps->ProjectionX("epem_4pi_Wolf_ps_px");
 epem_4pi_Wolf_ps_px->Draw();
 TH1* epem_4pi_Wolf_wt_ONLY_px = epem_4pi_Wolf_wt_ONLY->ProjectionX("epem_4pi_Wolf_wt_ONLY_px");
 epem_4pi_Wolf_wt_ONLY_px->Draw();

 
 TH1* epem_4pi_Wolf_wt_px1 = epem_4pi_Wolf_wt1->ProjectionX("epem_4pi_Wolf_wt_px");
 epem_4pi_Wolf_wt_px1->Draw();
 TH1* epem_4pi_Wolf_ps_px1 = epem_4pi_Wolf_ps1->ProjectionX("epem_4pi_Wolf_ps_px");
 epem_4pi_Wolf_ps_px1->Draw();

 TH1* epem_4pi_Wolf_wt_special_px = epem_4pi_Wolf_wt_special->ProjectionX("epem_4pi_Wolf_wt_special_px");
 epem_4pi_Wolf_wt_special_px->Draw();
 TH1* epem_4pi_Wolf_ps_special_px = epem_4pi_Wolf_ps_special->ProjectionX("epem_4pi_Wolf_ps_special_px");
 epem_4pi_Wolf_ps_special_px->Draw();
 
 //TH1* epem_4pi_Wolf_wt_DIV_px = epem_4pi_Wolf_wt_DIV_px->ProjectionX("epem_4pi_Wolf_wt_div_px");
 //epem_4pi_Wolf_wt_DIV_px->Draw();
 
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

 //hepemInvM->Scale(1./(2200000.));
 //hepemInvM->Scale(1./hepemInvM->GetEntries());   
 normalize(hepemInvM);
 hepemInvM->Draw();
 gr2->Draw("psame");

 TCanvas *can33=new TCanvas("brems_1D","brems_1D",1400,600);
 can33->Divide(3,3);

 can33->cd(1);
 hepemInvM->Draw();
 can33->cd(2);
 hydil->Draw();
 can33->cd(3);
 hptdil->Draw();

 can33->cd(4);
 hepemInvMCM->Draw();
 can33->cd(5);
 hydilCM->Draw();
 can33->cd(6);
 hptdilCM->Draw();
  
 
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
 hMypt_ps->SetFillColor(kCyan);
 hMypt_ps->Draw("ISO");
 can3->cd(11);
 hMypt->SetFillColor(kCyan);
 hMypt->Draw("ISO");
 can3->cd(12);
 hMypt_div->SetFillColor(kCyan);
 hMypt_div->Draw("ISO");



 //TFile *output = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/FINAL_CAL/Wolf_caln_PimP_Eng_Brem_FINAL_CHECK.root","RECREATE");
 //TFile *output = new TFile("Wolf_caln_PimP_Eng_Brem_FINAL_PS_CHECK.root","RECREATE");
 //TFile *output = new TFile("Wolf_caln_3D_HIST_METHOD_Eng_PimP_Brem_ANGULAR_Kine_Corr_FINAL_New_GiBUU_OCT2020.root","RECREATE"); //---------USE----THIS------PLEASE
 //TFile *output = new TFile("Wolf_caln_3D_HIST_METHOD_Eng_PimP_Brem_PHASESPACE_Kine_Corr_FINAL.root","RECREATE"); //---------USE----THIS------PLEASE
 //TFile *output = new TFile("Wolf_caln_3D_HIST_METHOD_Eng_PimN_PimN_Brem_PHASESPACE.root","RECREATE");
 //TFile *output = new TFile("Wolf_caln_PimN_Eng_Brem_FINAL_CS12_ANGULAR.root","RECREATE");
 TFile *output = new TFile("Wolf_caln_PimP_TEST_ANGULAR_dSigma_d3p.root","RECREATE");

 //TFile *output = new TFile("Wolf_caln_3D_HIST_METHOD_Eng_Pim_CARBON_Brem_PHASESPACE_3_3_Corr_FINAL.root","RECREATE");
 //TFile *output = new TFile("Caln_CARBON_TEST_RUN2.root","RECREATE");
  output->cd();

 hMypt_px->Write();
 hMypt_ps_px->Write(); 
 TH1F *Wolf_Brem = (TH1F*)hMypt_px->Clone();
 Wolf_Brem->Divide(hMypt_ps_px);
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

 epem_4pi_Wolf_wt_px1->Write();
 epem_4pi_Wolf_ps_px1->Write();
 TH1F *Wolf_Brem_3 = (TH1F*)epem_4pi_Wolf_wt_px1->Clone();
 Wolf_Brem_3->Divide(epem_4pi_Wolf_ps_px1);
 Wolf_Brem_3->SetName("Wolf_Brem3");
 Wolf_Brem_3->Draw("same,hist");
 Wolf_Brem_3->Write();

 epem_4pi_Wolf_wt_special_px->Write();
 epem_4pi_Wolf_ps_special_px->Write();
 TH1F *Wolf_Brem_special = (TH1F*)epem_4pi_Wolf_wt_special_px->Clone();
 Wolf_Brem_special->Divide(epem_4pi_Wolf_ps_special_px);
 Wolf_Brem_special->SetName("Wolf_Brem_special");
 Wolf_Brem_special->Draw("same,hist");
 Wolf_Brem_special->Write();


 epem_4pi_Wolf_wt_ONLY_px->Write();

 hMypt->Write();
 hMypt_ps->Write();
 PhS1->Write();
 PhS1_M->Write();
 PhS1_Ecm->Write();
 
 //epem_4pi_Wolf_wt->Write();       
 //epem_4pi_Wolf_ps->Write();       
/////////////////////////////////////////////////////////////////
 epem_4pi_Wolf_wt_DIV->Write();
 epem_4pi_Wolf_wt_DIV_px->SetName("Wolf_Brem_DIV");
 normalize(epem_4pi_Wolf_wt_DIV_px);
 epem_4pi_Wolf_wt_DIV_px->Write();
 /////////////////////////////////////////////////////////////////
 epem_4pi_Wolf_M_wt_DIV->Write();
 epem_4pi_Wolf_wt_special_DIV->Write();
 epem_4pi_Wolf_M_wt_DIV_px->SetName("Wolf_Brem_DIV_M");
 normalize(epem_4pi_Wolf_M_wt_DIV_px);
 epem_4pi_Wolf_M_wt_DIV_px->Write();
 /////////////////////////////////////////////////////////////////
 epem_4pi_Wolf_Ecm_wt_DIV->Write();
 epem_4pi_Wolf_Ecm_wt_DIV_px->SetName("Wolf_Brem_DIV_Ecm");
 normalize(epem_4pi_Wolf_Ecm_wt_DIV_px);
 epem_4pi_Wolf_Ecm_wt_DIV_px->Write();
 /////////////////////////////////////////////////////////////////

 
 //epem_4pi_Wolf_wt_DIV_px->Write();
 hMypt_div->Write();
 hMypt_div_px->SetName("Wolf_hMypt_div");
 normalize(hMypt_div_px);
 hMypt_div_px->Write(); //normalize(DS_invM);

 DS_inv->Write();
 DS_M->Write();
 DS_Ecm->Write();
 //DiffSec_saf->Write();
 
 output->Write();
 output->Close();
 /*
 TCanvas *h1Dvs3D=new TCanvas("h1Dvs3D","h1Dvs3D",1000,500);
 h1Dvs3D->cd();
 TH1D* epem_4pi_Wolf_wt_DIV_px = epem_4pi_Wolf_wt_DIV->ProjectionX("epem_4pi_Wolf_wt_DIV_px");
 epem_4pi_Wolf_wt_DIV_px->Draw("same,hist");
 Wolf_Brem_2->Draw("same,hist");
 */ 
  /*
   TCanvas *can=new TCanvas("brems_hist","brems_hist",1000,500);
   can->Divide(4,2);

   can->cd(1);
   gPad->SetLogy();   
   //hepemInvM->Scale(1./1000000.);   

   //normalize(hepemInvM);
   hepemInvM->Draw();
   //can->cd(2);
   //sigInvM->Scale(1./sigInvM->Integral());
   //sigInvM->Draw("colz"); 
   
   can->cd(2);
   hydil->Draw();
   
   can->cd(3);
   hptdil->Draw();

   can->cd(4);
   hEdil->Draw();


   can->cd(5);
   gPad->SetLogy();   
   //normalize(epemInvM);
   //epemInvM->Scale(1./epemInvM->Integral());   
   hepemInvM1->Draw();

   can->cd(6);
   hydil1->Draw();

   can->cd(7);
   hptdil1->Draw();

   can->cd(8);
   hEdil1->Draw();
   
   //can->cd(5);
   //yInvM->Draw("colz"); 
   
   //can->cd(6);
   //ptInvM->Draw("colz"); 
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



 
