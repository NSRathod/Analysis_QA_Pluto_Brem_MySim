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
#include <TLegend.h>
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
Double_t dsigdmpn_nopscorr_2M(float x);
Double_t dsigdmpn_nopscorr(float x);

void schafer_pp_test1(){
  ifstream f1a;
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/PN_Brem_Energy1GeV_1M.evt");
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PP_NEW_Bremsstrahlung.evt");
  //f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/TEST/Events_theoretical_Calc/Events_PPGamma_T1GeV.evt");
  f1a.open("/lustre/nyx/hades/user/nrathod/PLUTO/690/ANGULAR/pion690_BREM_ANG_22M.evt");
  

 const float PI= TMath::Pi();

 //double p_momentum = 1.6958;//GeV/c
 double p_momentum =0.690;//GeV/c
 //const float mpi=0.93827;
 const float mpi=0.13957;
 const float mp=0.93956;
 const float me=0.510998*0.001;

 double p_energy = sqrt(p_momentum*p_momentum + mp*mp);

 const float alpha = 1./137.;
 const float sqrts=1.5 ;
 //const float sigma_el=23.0;//mb 
 ///const float sqrts=2.32;
 const float sigma_el=23.0;//mb 
 const float A=alpha*alpha/(3.*PI*PI);
 const float A2=alpha*alpha/(4.*PI*PI);
 const float WA=alpha*alpha/(6.*PI*PI);
 const float WA2=alpha*alpha/(6.*PI*PI*PI);
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
   // hepemInvM->Sumw2();
   hepemInvM->SetLineColor(1);
   hepemInvM->SetMarkerColor(2);

   TH1F *hepemInvM_newcal = new TH1F("hepemInvM_newcal","hepemInvM_newcal",100,0.,0.6);
   // hepemInvM->Sumw2();
   hepemInvM_newcal->SetLineColor(2);
   hepemInvM_newcal->SetMarkerColor(2);
   
   TH1F *hepemInvM_haglin = new TH1F("hepemInvM_haglin","hepemInvM_haglin",100,0.,0.6);
   // hepemInvM->Sumw2();
   hepemInvM_haglin->SetLineColor(2);
   hepemInvM_haglin->SetMarkerColor(1);
   
   TH1F *hepemInvM_ps = new TH1F("hepemInvM_ps","hepemInvM_ps",100,0.,0.6);
   hepemInvM_ps->Sumw2();
   hepemInvM_ps->SetLineColor(3);
   hepemInvM_ps->SetMarkerColor(3);
   
   TH1F *hepemInvM_div = new TH1F("hepemInvM_div","hepemInvM_div",100,0.,0.6);
   hepemInvM_div->Sumw2();
   hepemInvM_div->SetLineColor(4);
   hepemInvM_div->SetMarkerColor(4);



   TH1F *hepemInvM_PLUTO_ph = new TH1F("hepemInvM_PLUTO_ph","hepemInvM_PLUTO_ph",111,0.0,0.555);
   TH1F *hepemInvM_PLUTO_wtph = new TH1F("hepemInvM_PLUTO_wtph","hepemInvM_PLUTO_wtph",111,0.0,0.555);
   TH1F *hepemInvM_PLUTO_wtph2 = new TH1F("hepemInvM_PLUTO_wtph2","hepemInvM_PLUTO_wtph2",111,0.0,0.555);
   TH1F *hepemInvM_PLUTO_wt = new TH1F("hepemInvM_PLUTO_wt","hepemInvM_PLUTO_wt",111,0.0,0.555);


   TH1F *hepemInvM_PLUTO1 = new TH1F("hepemInvM_PLUTO1","hepemInvM_PLUTO1",111,0.0,0.555);
   TH1F *hepemInvM_PLUTO1_ph = new TH1F("hepemInvM_PLUTO1_ph","hepemInvM_PLUTO1_ph",111,0.0,0.555);
   //TH1F *hepemInvM_PLUTO1_ph2 = new TH1F("hepemInvM_PLUTO_ph","hepemInvM_PLUTO_ph",111,0.0,0.555);

   TH1F *hepemInvM1 = new TH1F("hepemInvM_nocut","hepemInvM_nocut",111,0.0,0.555);

   TH1F *hepemMM = new TH1F("hepemMM","hepemMM",100,0.,1.2);
   TH1F *hOA = new TH1F("hOA","hOA",180, 0.0, 180);
   TH1F *hydil = new TH1F("hydil","hydil",100, -2.5, 2.5);
   TH1F *hptdil = new TH1F("hptdil","hptdil",100, 0, 0.5);
   TH1F *hydil1 = new TH1F("hydil_nocut","hydil_nocut",100, -2.5, 2.5);
   TH1F *hptdil1 = new TH1F("hptdil_nocut","hptdil_nocut",100, 0, 0.5);

   TH1F *hydil_ps = new TH1F("hydil_ps","hydil_ps",100, -2.5, 2.5);
   TH1F *hptdil_ps = new TH1F("hptdil_ps","hptdil_ps",100, 0, 0.5);
   TH1F *hydil_div = new TH1F("hydil_div","hydil_div",100, -2.5, 2.5);
   TH1F *hptdil_div = new TH1F("hptdil_div","hptdil_div",100, 0, 0.5);


   TH1F *hEdil = new TH1F("hEdil","hEdil",100, 0, 1.);
   TH1F *hEdil1 = new TH1F("hEdil_nocut","hEdil_nocut",100, 0, 1.);

   TH2F *sigInvM = new TH2F("sigInvM","sigInvM",100,0.,0.4,100,0,3);
   TH2F *yInvM = new TH2F("yInvM","yInvM",100,0.,0.4,100,-2.5,2.5);
   TH2F *ptInvM = new TH2F("ptInvM","ptInvM",100,0.,0.4,100,0,0.5);

   TH3F *epem_4pi_Wolf_wt = new TH3F("epem_4pi_Wolf_wt_Mypt","epem_4pi_Wolf_wt_Mypt", 111,0.0,0.555, 200,-3.0,4.0, 100, 0, 0.5);
   TH3F *epem_4pi_Wolf_wt_noph = new TH3F("epem_4pi_Wolf_wt_noph_Mypt","epem_4pi_Wolf_wt_noph_Mypt", 111,0.0,0.555, 200,-3.0,4.0, 100, 0, 0.5);
   TH3F *epem_4pi_Wolf_ps = new TH3F("epem_4pi_Wolf_ps_Mypt","epem_4pi_Wolf_ps_Mypt", 111,0.0,0.555, 200,-3.0,4.0, 100, 0, 0.5);
   TH3F *epem_4pi_Wolf_wt_ONLY = new TH3F("epem_4pi_Wolf_wt_Mypt_ONLY","epem_4pi_Wolf_wt_Mypt_ONLY", 200,0.0,1.0, 200,-3.0,4.0, 100, 0, 0.5);
   TH3F *epem_4pi_Wolf_wt_DIV = new TH3F("epem_4pi_Wolf_wt_Mypt_DIV","epem_4pi_Wolf_wt_Mypt_DIV", 111,0.0,0.555, 200,-3.0,4.0, 100, 0, 0.5);   
   TH3F *epem_4pi_Wolf_wt_noph_DIV = new TH3F("epem_4pi_Wolf_wt_noph_Mypt_DIV","epem_4pi_Wolf_wt_noph_Mypt_DIV", 111,0.0,0.555, 200,-3.0,4.0, 100, 0, 0.5);   
   TH1F *DiffSec = new TH1F("DiffSec","DiffSec",200,0.,1.0); DiffSec->Sumw2();
   TH1F *DiffSec_saf = new TH1F("DiffSec_saf","DiffSec_saf",200,0.,1.0); DiffSec_saf->Sumw2();
   TH1F *DiffSec_wolf = new TH1F("DiffSec_wolf","DiffSec_wolf",200,0.,1.0); DiffSec_wolf->Sumw2();
   
   TF1*  Brempn_nops=new TF1 ("Brempn_nops","dsigdmpn_nopscorr(x)",0.005,0.555);   
   TF1*  Brempn_nops_2M=new TF1 ("Brempn_nops","dsigdmpn_nopscorr_2M(x)",0.005,0.555);   

   float sig[8]={0.};
   float nb[8]={0.};
   TGraph *gr1;
   //float nbEv=1000000.;
   //float nbEv=595609.;
   float s=sqrts*sqrts;


   
   //TFile *f1 = new TFile("Gibuu_Calc.root","RECREATE");    
   //ifstream myfile ("TEST.txt");
   float ii = 0;
   ofstream f1;
   f1.open ("Schafers_approach_pp_Brem.evt");
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
  
  lv->SetPxPyPzE(c2,c3,c4,c1);	 
  lv1->SetPxPyPzE(d2,d3,d4,d1);	 
  flag1=0;
  flag2=0;
  
  if(c9==2){
    p3.SetXYZ(lv->Px(),lv->Py(),lv->Pz());
    ep->SetVectM(p3, me);
    flag1=1;//
     }		      

  
  if(d9==3){
    p4.SetXYZ(lv1->Px(),lv1->Py(),lv1->Pz());
    em->SetVectM(p4, me);
    flag2=1;	      
      }
  

    *mm=*beam-*ep;
    *invm = *ep + *em;
    
    float Invmass = invm->M();
    float InvEn = invm->E();
    //float Ep=ep->E();
    //float Em=em->E();
    //float en1=InvEn*sqrt(InvEn*InvEn-Invmass*Invmass);
    //cout<<InvEn<<" "<<Ep<<" "<<Em<<" "<<Ep+Em<<endl;
    //cout<<"Mass :: "<<Invmass<<endl;
    float InvMom = invm->P();
    float pt = invm->Pt();
    float pt2 = invm->Pt()*invm->Pt();
    float y = invm->Rapidity();
    float y_ep = ep->Rapidity();
    float y_em = em->Rapidity();
    float TH = invm->Theta() * TMath::RadToDeg();
    float theta = TMath::Sin(TH);
    float phi = invm->Phi();
    float dOmega = 2*PI*theta*TH*phi;
    //cout<<InvEn<<" "<<InvEnCM<<endl;
    //cout<<"theta: "<<theta<<"   "<<"TH: "<<TH<<"   "<<"phi :"<<phi<<"   "<<"Omega : "<<dOmega<<endl;
    
    float pt_ep = ep->Pt();
    float pt_em = em->Pt();
    //float missmass = mm->M();

    float momep=ep->P();
    float momem=em->P();
    float m2=(momep+momem)*(momep+momem);
    
    ep->Boost(-(*beam).BoostVector() );
    em->Boost(-(*beam).BoostVector() );
    *invmCM = *ep + *em;
    
    float InvMomCM = invmCM->P();
    float InvEnCM = invmCM->E();
    float yCM = invmCM->Rapidity();
    float ptCM = invmCM->Pt();

    float momepCM=ep->P();
    float momemCM=em->P();
    float m2CM=(momepCM+momemCM)*(momepCM+momemCM);
/////////----------Addition--------------------------------------------
    //--------------ALL required phasespace scenarios -----------------
    float Energy = sqrt(pow(Invmass,2)+pow(pt,2))*TMath::CosH(y);
    float s2=s+(Invmass*Invmass)-(2.*Energy*sqrts);
    float s2_M=s+(Invmass*Invmass)-(2.*Invmass*sqrts);
    float s2_Ecm=s+(Invmass*Invmass)-(2.*InvEnCM*sqrts);
    
    float R2s=sqrt(1.-(((mp+mpi)*(mp+mpi))/s));
    float R2s2=sqrt(1.- (((mp+mpi)*(mp+mpi))/s2)); 
    float R2s2_M=sqrt(1.- (((mp+mpi)*(mp+mpi))/s2_M)); 
    float R2s2_Ecm=sqrt(1.- (((mp+mpi)*(mp+mpi))/s2_Ecm)); 

    //float R2s2_M= (2*Qcm_s2_M)/sqrt(s2_M);
    //float R2s2_Ecm= (2*Qcm_s2_Ecm)/sqrt(s2_Ecm);
/////////----------Addition---------------------------------------------

    float sigma_quasiel_PLUTO= 2*((s-((mpi+mp)*(mpi+mp)))/((mpi+mp)*(mpi+mp)))*18.0;    //SCHAFER's effective SIGMA
    float sigma_quasiel_PLUTO1= ((s-((mpi+mp)*(mpi+mp)))/(2*(mpi*mpi)))*18.0;           //WOLF's effective SIGMA
    double RATIO_sig_bar = sigma_quasiel_PLUTO/sigma_quasiel_PLUTO1;  // RATIO OD ABOVE 2 effective sigma
    //float sigma_quasiel_PLUTO= 2.*sigma_el*((s/(4.*mp*mp)) -1.);
    //float dsigma=A*(sigma_quasiel/(Invmass*InvEn*InvEn))*(R2s2/R2s)*y*pt*pt;
    float dsigma_PLUTO=A*(sigma_quasiel_PLUTO1/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass)*2*Invmass;   //-----SCHAFER's FORMULA (dSigma/dM)
    float dsigma_PLUTO_ph2=A*(sigma_quasiel_PLUTO1/(Invmass*Invmass))*log((Invmass+0.35)/Invmass)*2*Invmass*(R2s2_M/R2s);
    float dsigma_PLUTO_ph=A*(sigma_quasiel_PLUTO1/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass)*2*Invmass*(R2s2_M/R2s);
    //SCHAFER's FORMULA (dSigma/dM)+PHASESPACE
    float dsigma_PLUTO1=A*(sigma_quasiel_PLUTO1/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass);  ////-----SCHAFER's FORMULA (dSigma/dM2)
  float dsigma_PLUTO_Method2=WA2*(sigma_quasiel_PLUTO1/(Invmass*pow(InvEnCM,3)))*(R2s2/R2s);  ////-----SCHAFER's FORMULA WORK_IN_PROGRESS (dSigma/d3p)


    float dsigmaW=WA*(sigma_quasiel_PLUTO1/(Invmass*Energy*Energy))*(R2s2_M/R2s);  /////-----WOLF's FORMULA + noPhasespace
    float dsigmaW1=WA*(sigma_quasiel_PLUTO1/(Invmass*Energy*Energy));//*(1/(2*Invmass));      /////-----WOLF's FORMULA + noPhasespace
    float dsigmaW2=WA2*(sigma_quasiel_PLUTO1/(Invmass*Energy*Energy))*InvMom*(R2s2/R2s);   //-----WOLF's FORMULA (dSigma/d3p)//============WORK--IN--PROGRESS
    //cout <<"P.Schafer sig = "<< sigma_quasiel_PLUTO<<"   "<<"P.W0lf's sig = "<<sigma_quasiel_PLUTO1 << "      "<<RATIO_sig_bar<<endl;
    //Double_t Opang = openingangle(*ep,*em);
    //Double_t Opang_deg = Opang*TMath::RadToDeg();

    //if (TMath::IsNaN(dsigmaW))continue;
    std::cout.precision(20);
    //ii = ii+1; 
    //cout<<":::::"<<ii<<endl;
    //cout<<InvEn<<" "<<InvEnCM<<endl;
    //cout<<"y: "<<y<<" "<<yCM<<endl;
    //cout<<"pt: "<<pt<<" "<<ptCM<<endl;
    //cout<<"dSigma="<<dsigdmpn_nopscorr(Invmass)<<endl;
    
    //hepemInvM_PLUTO_ph->Fill(Invmass);

    if (TMath::IsNaN(dsigma_PLUTO_ph))continue;
    hepemInvM_PLUTO_wtph->Fill(Invmass,dsigma_PLUTO_ph);
    hepemInvM_PLUTO_wtph2->Fill(Invmass,dsigma_PLUTO_ph2);
    hepemInvM_PLUTO_wt->Fill(Invmass,dsigma_PLUTO);
    hepemInvM_PLUTO_ph->Fill(Invmass);
    hepemInvM_PLUTO1->Fill(Invmass,dsigma_PLUTO1);
    hepemInvM_PLUTO1_ph->Fill(Invmass);

    if (TMath::IsNaN(dsigmaW)) continue;
    ii = ii+1; 
    cout<<":::::"<<ii<<endl;
    epem_4pi_Wolf_wt->Fill(Invmass,y,pt2,dsigmaW);
    epem_4pi_Wolf_wt_noph->Fill(Invmass,y,pt2,dsigmaW1);
    epem_4pi_Wolf_ps->Fill(Invmass,y,pt2);

    /*    if (TMath::IsNaN(dsigmaW2)) continue;
    epem_4pi_Wolf_wt_2->Fill(Invmass,y,pt,dsigmaW2);
    epem_4pi_Wolf_wt_noph->Fill(Invmass,y,pt,dsigmaW1);
    epem_4pi_Wolf_ps->Fill(Invmass,y,pt);
    */ //----TO-DO--LIST

    
    //DiffSec->Fill(dsigmaW);
    //DiffSec_saf->Fill(dsigma_PLUTO1);
    //DiffSec_wolf->Fill(dsigma_PLUTO);
    //DS_invM->Fill(dsigma_PLUTO);
    
    /*    
   cout<<p1a<<" "<<wa<<" "<<w1a<<" "<<w2a<<" "<<w3a<<endl;
   cout<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a13*dsigma_PLUTO1<<endl;
   cout<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<" "<<b11<<" "<<b12<<" "<<b13*dsigma_PLUTO1<<endl;
   cout<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<" "<<c11<<" "<<c12<<" "<<c13*dsigma_PLUTO1<<endl;
   cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<" "<<d7<<" "<<d8<<" "<<d9<<" "<<d10<<" "<<d11<<" "<<d12<<" "<<d13*dsigma_PLUTO1<<endl;
    
    
   f1<<p1a<<" "<<wa<<" "<<w1a<<" "<<w2a<<" "<<w3a<<endl;
   f1<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<" "<<a7<<" "<<a8<<" "<<a9<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a13*dsigma_PLUTO1<<"\n";
   f1<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" "<<b10<<" "<<b11<<" "<<b12<<" "<<b13*dsigma_PLUTO1<<"\n";
   f1<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c6<<" "<<c7<<" "<<c8<<" "<<c9<<" "<<c10<<" "<<c11<<" "<<c12<<" "<<c13*dsigma_PLUTO1<<"\n";
   f1<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<" "<<d7<<" "<<d8<<" "<<d9<<" "<<d10<<" "<<d11<<" "<<d12<<" "<<d13*dsigma_PLUTO1<<"\n";
    */

    //float s2=s+(Invmass*Invmass)-(2.*InvEnCM*sqrts);
    //float s2=s+(Invmass*Invmass)-(2.*en1*sqrts);
   }
f1.close();


 TH1* epem_4pi_Wolf_wt_px = epem_4pi_Wolf_wt->ProjectionX("epem_4pi_Wolf_wt_px");
 epem_4pi_Wolf_wt_px->Draw();
 TH1* epem_4pi_Wolf_ps_px = epem_4pi_Wolf_ps->ProjectionX("epem_4pi_Wolf_ps_px");
 epem_4pi_Wolf_ps_px->Draw();

//----------------------------------------------------------------------------------------------
 epem_4pi_Wolf_wt->Scale(1/epem_4pi_Wolf_wt->Integral()); 
 epem_4pi_Wolf_wt_noph->Scale(1/epem_4pi_Wolf_wt_noph->Integral()); 
   epem_4pi_Wolf_ps->Scale(1/epem_4pi_Wolf_ps->Integral()); 

   epem_4pi_Wolf_wt_DIV->Divide(epem_4pi_Wolf_wt,epem_4pi_Wolf_ps,1.,1.);
 epem_4pi_Wolf_wt_DIV->Draw("ISO");

 epem_4pi_Wolf_wt_noph_DIV->Divide(epem_4pi_Wolf_wt_noph,epem_4pi_Wolf_ps,1.,1.);
 epem_4pi_Wolf_wt_noph_DIV->Draw("ISO");
//----------------------------------------------------------------------------------------------
 TH3F *Wolf_Brem_DIV_3D = (TH3F*)epem_4pi_Wolf_wt->Clone();
 Wolf_Brem_DIV_3D->Divide(epem_4pi_Wolf_ps); 
//----------------------------------------------------------------------------------------------
 TH1* Wolf_Brem_DIV_3D_px = Wolf_Brem_DIV_3D->ProjectionX("Wolf_Brem_DIV_3D_px");                     //----Normal_Division
 TH1* epem_4pi_Wolf_wt_DIV_px = epem_4pi_Wolf_wt_DIV->ProjectionX("epem_4pi_Wolf_wt_DIV_px");         //----Division
 TH1* epem_4pi_Wolf_wt_noph_DIV_px = epem_4pi_Wolf_wt_noph_DIV->ProjectionX("epem_4pi_Wolf_wt_noph_DIV_px");         //----Division


 /*
 while (!f1a.eof()){
   // cout<<a1<<endl;                                      
 }*/
    //float R2s=sqrt(1.-((mp+mp)*(mp+mp)/s));
    //float R2s2=sqrt(1.- ( (mp+mp)*(mp+mp)/s2) );
    //float R=sqrt((s*(s2-4*mp*mp))/(s2*(s-4*mp*mp)));
    //float sigma_quasiel= ((s-(mp+mp)*(mp+mp))/(2.*mp*mp))*sigma_el;
/*
   for (int i=0;i<50;i++){

     float mass=0.01+i*0.01;



    float sigma_quasiel= ((sqrts*sqrts-((mpi+mp)*(mpi+mp)))/(2.*mpi*mpi))*sigma_el;
    //float sigma_quasiel_PLUTO= 2.*sigma_el*((s/(4.*mp*mp)) -1.);
    //float dsigma=A*(sigma_quasiel/(Invmass*InvEn*InvEn))*(R2s2/R2s)*y*pt*pt;
    float dsigma=A*(sigma_quasiel/(mass*mass))*log((sqrts-(mpi+mp))/mass);
    float dsigma_new=A*(sigma_quasiel/(mass*mass))*log((sqrts-(mpi+mp))/mass)*2*mass;

     float sigma_quasiel= 2.*sigma_el*((s/(4.*mp*mp)) -1.);
    //float sigma_quasiel1= 2.*sigma_el*((s/(4*pow((mp+mpi),2))) -1.);
    float sigma_quasiel1= sigma_el*((s - ((mpi + mp)*(mpi + mp)))/(2*(mpi*mpi)));
    //float dsigma=A*(sigma_quasiel/(Invmass*InvEn*InvEn))*(R2s2/R2s)*y*pt*pt;
    float dsigma=A*(sigma_quasiel1/(mass*mass))*log((sqrts-(mpi+mp))/mass);
    //float dsigma=A*(sigma_quasiel/(mass*mass))*log((sqrts-2.*mp)/mass);
    float dsigma_new=A*(sigma_quasiel1/(mass*mass))*log((sqrts-(mp+mpi))/mass)*2*mass;
    

    float sigma_quasielastic= sigma_el*((s - ((mpi + mp)*(mpi + mp)))/(2*(mpi*mpi)));
    float L1 = log(((sqrts - (mp+mpi)) + sqrt(pow((sqrts - (mp+mpi)),2)-pow(mass,2))) / ((sqrts - (mp+mpi)) - sqrt(pow((sqrts - (mp+mpi)),2)-pow(mass,2))));
    float L2 = 2*sqrt(pow((sqrts - (mp+mpi)),2)-pow(mass,2))/(sqrts - (mp+mpi)));
    float dsigmaDM2=A2*(L1 - L2)*sigma_quasielastic*(1 / pow(mass,2))*2*mass;

    if (TMath::IsNaN(dsigmaDM2))continue;

    //dsigma=dsigma/(Invmass*Invmass);

    cout << dsigma <<"  "<<dsigma_new<<endl;
    //cout<<"s2: "<<s2<<" R2s2: "<<R2s2<<" R2s: "<<R2s<<"  "<<sigma_quasiel<<"  "<<endl;
    //cout<<"::: "<<R2s2<<"  "<<((mpi+mp)*(mpi+mp)/s2)<<endl;

    //if (((mpi+mp)*(mpi+mp)/s2)>1) cout<<Invmass<<" "<<InvEn<<endl;
    //cout<<Invmass<<" "<<dsigma<<"  "<<s2<<"  "<<R2s2<<"  "<<R2s<<"  "<<sigma_quasiel<<"  "<<s2<<endl;
    //cout<<Invmass<<" "<<dsigma<<"  "<<y<<"  "<<pt<<endl;

    //hEdil1->Fill(InvEn);
    //hydil1->Fill(y);
    //hptdil1->Fill(pt);
    hepemInvM1->Fill(mass);

    
    //if (((mp+mp)*(mp+mp)/s2)>1)continue;
    //  if (TMath::IsNaN(dsigma))continue;
    
    //sigInvM->Fill(Invmass,dsigma);
    //yInvM->Fill(Invmass,y);
    //ptInvM->Fill(Invmass,pt);

    //hepemInvM->Fill(Invmass,dsigma);
    //hydil->Fill(y);
    //hptdil->Fill(pt);
    //hEdil->Fill(InvEn);

    hepemInvM->Fill(mass,dsigma);
    hepemInvM_newcal->Fill(mass,dsigma_new);
    //hepemInvM_haglin->Fill(mass,dsigmaDM2);
    //hydil->Fill(y,dsigma/nbEv);
    //hptdil->Fill(pt,dsigma/nbEv);

    
    //hepemInvM_ps->Fill(Invmass,1./nbEv);
    //hydil_ps->Fill(y,1./nbEv);
    //hptdil_ps->Fill(pt,1./nbEv);

    
  
    }*/

/*
 for (int i=0;i<50;i++){

   float mass=0.01+i*0.01;
   float sigma_quasiel= ((sqrts*sqrts-((mpi+mp)*(mpi+mp)))/(2.*mpi*mpi))*sigma_el;
   //float sigma_quasiel= 2.*sigma_el*((s/(4.*mpi*mpi)) -1.);
   //float dsigma=A*(sigma_quasiel/(Invmass*InvEn*InvEn))*(R2s2/R2s)*y*pt*pt;
   //float dsigma=A*(sigma_quasiel/(mass*mass))*log((sqrts-2.*mpi)/mass);
   //float dsigma=A*(sigma_quasiel/m2CM)*log((sqrts-2.*mp)/Invmass);
   float dsigma11=A*(sigma_quasiel_PLUTO/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass)*2*Invmass;

   hepemInvM1->Fill(mass);
   
   float sigma_quasielastic= ((sqrts*sqrts-((mpi+mp)*(mpi+mp)))/(2.*mpi*mpi))*sigma_el;
   float T1 = log( (sqrts - (mp+mpi)) + sqrt( pow((sqrts - (mp+mpi)),2) - pow((mass),2)  ));
   float T2 = log( (sqrts - (mp+mpi)) - sqrt( pow((sqrts - (mp+mpi)),2) - pow((mass),2)  ));
   float T3 = (2*sqrt( pow((sqrts - (mp+mpi)),2) - pow((mass),2))) / (sqrts - (mp+mpi));

   float dsigmaDM2=A2*((T1 - T2) - T3)*sigma_quasielastic*(1 / mass)*2;     

   cout << sigma_quasielastic<<"  "<<T1<<"  "<<T2<<"  "<<T3<"  "<<dsigmaDM2<<endl;
   //if (((mp+mp)*(mp+mp)/s2)>1)continue;
     //if (TMath::IsNaN(dsigmaDM2))continue;

hepemInvM->Fill(mass,dsigma);
//hepemInvM_newcal->Fill(mass,dsigma_new); 
//hepemInvM_haglin->Fill(mass,dsigmaDM2); 

 }
*/      
 //normalize(hepemInvM_PLUTO_ph);
   //normalize(hepemInvM_PLUTO1);
   //normalize(hepemInvM_PLUTO);
   //normalize(mass); 
   //hepemInvM_div->Divide(hepemInvM,hepemInvM_ps,1.,1.);
   //hydil_div->Divide(hydil,hydil_ps,1.,1.);
   //hptdil_div->Divide(hptdil,hptdil_ps,1.,1.);

 //  TGraph *gr=new TGraph(8,invM,sigfin);
 // gr->SetMarkerStyle(21);
 // gr->SetMarkerSize(1.);
 
 //TCanvas *can1=new TCanvas("brems_theory1","brems_theory1");
 //can1->cd();
 //gPad->SetLogy();
 //gr->Draw("al");

   float xSchaMathe2[]={0.052089, 0.0277064, 0.0157588, 0.0102732, 0.00780263, 0.00562184, 0.00417451, 0.00322255, 0.00238727, 0.00193423, 0.00146358, 0.00118192, 0.000835694, 0.000703856, 0.000588379, 0.000488188, 0.000415006, 0.000334985, 0.000276522, 0.000233062, 0.000193997, 0.000168262, 0.000137592, 0.00010779, 0.000085356, 0.0000673488, 0.0000543745, 0.0000415802, 0.0000311007 };
   
   //float xSchaMathe[]={0.0430301, 0.0228879, 0.0130181, 0.00848656, 0.00644565, 0.00464413, 0.00344851, 0.0026621, 0.00197209, 0.00159784, 0.00120904, 0.000976368, 0.000690356, 0.000581446, 0.000486052, 0.000403286, 0.000342831, 0.000276726, 0.000228431, 0.00019253, 0.000160258, 0.000138999, 0.000113663, 0.0000890442, 0.0000705115, 0.000055636, 0.0000449181, 0.0000343488, 0.0000256919};
 
 float xScha[]={0.04427643,0.057239614,0.07163615,0.08459933,0.09396338,0.106249265,0.118519396,0.13008073, 0.14454027,0.1553613,0.17054538,0.18281552,0.20378815,0.21461704, 0.22618626,0.23848002, 0.24933255,0.2637921,0.27679464, 0.28835598, 0.30064973,0.31005317,0.32306358,0.3382477,0.35198268,0.3649931,0.37582988,0.38809213,0.3996456 };

 float yScha[]={0.03325728,0.020907206,0.011987566,0.0075359894,0.005444284,0.004499853,0.00339556,0.0028079145,0.002317377,0.0017504091,0.0014439011,0.0010895585,8.553756E-4,	6.761951E-4,5.852167E-4,5.062283E-4,4.587523E-4,3.786091E-4,2.9885562E-4,2.471348E-4,2.1377832E-4,1.9392121E-4,1.6020199E-4,1.3214958E-4,1.09117274E-4,	9.0143854E-5,7.458019E-5,5.3773012E-5,4.2487834E-5};

 //---------------REAL POINTS FROMWOLF PAPER---------------------------------//
 float xWolf[]={0.05565393,0.06698674,0.07376166,0.08132172,0.089637786,0.09870986,0.10778921,0.11763914,0.12673306, 0.1350564, 0.14492819, 0.15325882, 0.16236003, 0.17146122, 0.18131845, 0.18966363,0.19799428,0.20861477, 0.21847199, 0.2275659,0.23819369, 0.2480509, 0.25944197, 0.2692846,0.28066114,0.29126707,0.30188027,0.31171563,0.32003897,0.32987434, 0.3389537,0.3502938, 0.3578247,0.3653848, 0.37217426,0.3804612,0.38722154,0.3939819,0.39848882};

 float yWolf[]={0.002029871, 0.0015116982, 0.0011904291, 9.900433E-4,8.0835074E-4, 6.4794865E-4, 5.289451E-4, 4.396775E-4,3.7227408E-4, 3.0955495E-4, 2.7179986E-4, 2.301725E-4, 1.9847752E-4, 1.7114694E-4, 1.4488454E-4,  1.2725792E-4,1.0776782E-4, 9.121501E-5, 7.721812E-5, 6.538044E-5, 5.635784E-5, 4.770976E-5, 4.111858E-5, 3.3560813E-5, 2.7887192E-5, 2.2757444E-5,1.891348E-5, 1.5157812E-5, 1.26040895E-5, 1.0101283E-5, 8.246061E-6, 6.2542103E-6,  4.835106E-6, 4.0212085E-6, 3.2843816E-6, 2.49277E-6,  1.8926141E-6, 1.436951E-6, 1.195901E-6};



 TGraph *gr22 = new TGraph(36,xWolf,yWolf);
 gr22->SetMarkerStyle(21);
 gr22->SetMarkerSize(1.);
 gr22->SetMarkerColor(1);   
 gr22->SetName("Wolf_points");
 TCanvas *can4_brem=new TCanvas("brems_Wold_DIV-hist","brems_Wold_DIV-hist");
 can4_brem->cd();
 epem_4pi_Wolf_wt_DIV_px->Draw("same");  
 gr22->Draw("P,same");

 TGraph *gr2 = new TGraph(29,xScha,yScha);
 gr2->SetMarkerStyle(21);
 gr2->SetMarkerSize(1.);
 gr2->SetMarkerColor(2);
 gr2->SetName("Schafer_points");
 
 TGraph *gr3 = new TGraph(29,xScha,xSchaMathe2);
 gr3->SetMarkerStyle(20);
 gr3->SetMarkerSize(1.);
 gr3->SetMarkerColor(3);
 gr3->SetName("Saf_Mathe_points");

 
 TLegend *leg = new TLegend(0.17,0.65,0.48,0.88);
 //leg->SetHeader("The Legend Title","C");
 leg->AddEntry(hepemInvM,"schafer","l");
 leg->AddEntry(hepemInvM_div,"schafer/Pluto","l");
 leg->AddEntry(hepemInvM_ps,"Pluto ps","l");
 leg->AddEntry(gr2,"schafer paper@1GeV","p");

 leg->SetTextFont(132);
 leg->SetTextSize(0.05);
 leg->SetBorderSize(0);

 //TCanvas *can2=new TCanvas("brems_th_wolf_en1","brems_th_wolf_en1");
 TCanvas *can2=new TCanvas("brems_th_Schafer","brems_th_Schafer");

 can2->cd();
 gPad->SetLogy();
 //hepemInvM->Scale(1./(1000000.*6));
 // hepemInvM->Scale(1./(1000000.));

 //hepemInvM->Scale(1./(5.));
 //hepemInvM->Scale(1./hepemInvM->GetEntries());   

 //normalize(hepemInvM);
 //normalize(hepemInvM_div);
 //normalize(hepemInvM_ps);

 //hepemInvM_ps->Draw();
 //hepemInvM_div->Draw("same");
 hepemInvM->Draw("same");
 gr2->Draw("psame");
 gr3->Draw("psame");


 TCanvas *can22=new TCanvas("brems_th_Schafer_vs_wolf","brems_th_Schafer_vs_wolf");

 can22->cd();
 gPad->SetLogy();
 //hepemInvM_newcal->SetLineColor(2);
 hepemInvM_newcal->Draw("same");
 gr22->Draw("psame");
hepemInvM_haglin->Draw("same");
//leg->Draw();

 TCanvas *can3=new TCanvas("brems_EpEm_dsigmadM_Schafer","brems_EpEm_dsigmadM_Schafer");
 can3->cd();
 gPad->SetLogy();
 hepemInvM_PLUTO1->Divide(hepemInvM_PLUTO1_ph);
 hepemInvM_PLUTO1->Draw("hist");
 gr2->Draw("P,same");

 TCanvas *can4=new TCanvas("brems_EpEm_Schafer","brems_EpEm_Schafer");
 can4->cd();
 gPad->SetLogy();
 hepemInvM1->Draw("hist");
 /*
 TCanvas *can5=new TCanvas("brems_EpEm_dsigmadM2_Schafer","brems_EpEm_dsigmadM2_Schafer");
 can5->cd();
 gPad->SetLogy();
 hepemInvM_PLUTO1->Divide(hepemInvM_PLUTO1_ph);
 hepemInvM_PLUTO1->Draw("hist");
 gr2->Draw("psame");
 */
 
 /*
 TCanvas *can3=new TCanvas("brems_dsig","brems_dsig");
  can3->cd();
  gPad->SetLogy();
   gPad->SetLogz();
  //sigInvM->Scale(1./100000.);   
  //hepemInvM->Scale(1./hepemInvM->Integral());   
  sigInvM->Draw("colz");
 */
 TCanvas *can51=new TCanvas("brems_EpEm_dsigmadM_Schafer_ph","brems_EpEm_dsigmadM_Schafer_ph");
 can51->cd();
 gPad->SetLogy();
 hepemInvM_PLUTO_ph->Draw("hist");
 
 /*
 TCanvas *can3=new TCanvas("brems_dsig","brems_dsig");
  can3->cd();
  gPad->SetLogy();
   gPad->SetLogz();
  //sigInvM->Scale(1./100000.);   
  //hepemInvM->Scale(1./hepemInvM->Integral());   
  sigInvM->Draw("colz");
 */
 
 TCanvas *can25=new TCanvas("brems_EpEm_dsigmadM_haglin","brems_EpEm_dsigmadM_haglin");
 can25->cd();
 gPad->SetLogy();
 hepemInvM_haglin->Draw("hist");
 

 //TFile *output = new TFile("Schafer_calculation_PP_1GeV_Bremsstrahlung_FINAL_Const_PI2_CS20_M2.root","RECREATE");
 //TFile *output = new TFile("Schafer_calculation_PP_1GeV_Bremsstrahlung_FINAL_OCT2020.root","RECREATE");
 //TFile *output = new TFile("Schafer_Wolf_PP_1GeV_TEST_RUN.root","RECREATE"); 
 TFile *output = new TFile("Schafer_Wolf_PIMP_1GeV_TEST_RUN_18mb_mul2M_WOLF.root","RECREATE");
 //TFile *output = new TFile("Schafer_Wolf_PN_1GeV_TEST_RUN_18mb_nomass_in_wtgs.root","RECREATE");
 output->cd();
 hepemInvM_PLUTO_ph->Write();
 //hepemInvM_PLUTO->Write();
 hepemInvM_PLUTO1->Write();
 hepemInvM1->Write();
 hepemInvM_haglin->Write();
 gr2->Write();
 gr22->Write();


 epem_4pi_Wolf_wt_px->Write();
 epem_4pi_Wolf_ps_px->Write();
 TH1F *Wolf_Brem_2 = (TH1F*)epem_4pi_Wolf_wt_px->Clone();
 Wolf_Brem_2->Divide(epem_4pi_Wolf_ps_px);
 Wolf_Brem_2->SetName("Wolf_Brem2");
 Wolf_Brem_2->Draw("same,hist");
 Wolf_Brem_2->Write();
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
 epem_4pi_Wolf_wt_DIV->Write();
 epem_4pi_Wolf_wt_DIV_px->SetName("Wolf_Brem_DIV");  
 //normalize(epem_4pi_Wolf_wt_DIV_px);
 epem_4pi_Wolf_wt_DIV_px->Write();  

 epem_4pi_Wolf_wt_noph_DIV->Write();
 epem_4pi_Wolf_wt_noph_DIV_px->SetName("Wolf_Brem_noph_DIV");  
 //normalize(epem_4pi_Wolf_wt_noph_DIV_px);
 epem_4pi_Wolf_wt_noph_DIV_px->Write();
 
 //normalize(DS_invM);
 DiffSec->Write();
 DiffSec_saf->Write();
 DiffSec_wolf->Write();
  
 Wolf_Brem_DIV_3D->SetName("Wolf_Brem_DIV_3D_new");
 Wolf_Brem_DIV_3D->Write();
 normalize(Wolf_Brem_DIV_3D_px);
 Wolf_Brem_DIV_3D_px->Write();


 normalize(hepemInvM_PLUTO_ph);
 normalize(hepemInvM_PLUTO_wtph);
 normalize(hepemInvM_PLUTO_wtph2);
 normalize(hepemInvM_PLUTO_wt);

 hepemInvM_PLUTO_ph->Write();
 hepemInvM_PLUTO_wtph->Divide(hepemInvM_PLUTO_ph);
 hepemInvM_PLUTO_wtph2->Divide(hepemInvM_PLUTO_ph);
 //normalize(hepemInvM_PLUTO_wtph);
 //hepemInvM_PLUTO_wtph->Scale(0.001);
 hepemInvM_PLUTO_wtph->Write(); 
 hepemInvM_PLUTO_wtph2->Write(); 
 hepemInvM_PLUTO_wt->Divide(hepemInvM_PLUTO_ph);
 //normalize(hepemInvM_PLUTO_wt);
 //hepemInvM_PLUTO_wt->Scale(0.009);
 hepemInvM_PLUTO_wt->Write();

 TCanvas *can_new=new TCanvas("brems_EpEm_Schafer_NEWEST","brems_EpEm_Schafer_NEWEST");
 can_new->cd();
 gPad->SetLogy();
 //hepemInvM_PLUTO_wtph->Scale(0.000001);
 //hepemInvM_PLUTO_wtph->Draw("hist");
 //hepemInvM_PLUTO_wt->Draw("hist,same");
 //Brempn_nops->Draw("same");
 //Brempn_nops_2M->Draw("same");

 epem_4pi_Wolf_wt_noph_DIV_px->Draw("hist,same");
 epem_4pi_Wolf_wt_DIV_px->Draw("hist,same");
 gr22->Draw("P,same");
 //hepemInvM_PLUTO1->Draw("hist,same");
 output->Write();
 output->Close();


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



 
Double_t dsigdmpn_nopscorr_2M(float x){

 const float PI1= TMath::Pi();
 //const float mpi1=0.93827;
 const float mpi1=0.13957;
 const float mp1=0.93956;

 const float alpha1 = 1./137.;
 const float sqrts1=1.5;//2.32;
 const float s1 = sqrts1*sqrts1;
 
 const float A1=(alpha1*alpha1)/(3.*PI1*PI1);
 double_t facsig=2*18*((s1-((mpi1+mp1)*(mpi1+mp1)))/((mpi1+mp1)*(mpi1+mp1)));
 //cout << "facsig: " << facsig <<"   "<<"x : "<< x <<endl;
 return A1 * (facsig / (x*x)) * log((sqrts1-(mpi1+mp1))/x)*2*x;

}
Double_t dsigdmpn_nopscorr(float x){

 const float PI2= TMath::Pi();
 //const float mpi2=0.93827;
 const float mpi2=0.13957;
 const float mp2=0.93956;

 const float alpha2 = 1./137.;
 const float sqrts2=1.5;//2.32;
 const float s2 = sqrts2*sqrts2;
 
 const float A2=(alpha2*alpha2)/(3.*PI2*PI2);
 double_t facsig=2*18*((s2-((mpi2+mp2)*(mpi2+mp2)))/((mpi2+mp2)*(mpi2+mp2)));
 //cout << "facsig: " << facsig <<"   "<<"x : "<< x <<endl;
 return A2 * (facsig / (x*x)) * log((sqrts2-(mpi2+mp2))/x);

}
/*
float sigma_quasiel_PLUTO= 2*((s-((mpi+mp)*(mpi+mp)))/((mpi+mp)*(mpi+mp)))*18.0;    //SCHAFER's effective SIGMA
float sigma_quasiel_PLUTO1= ((s-((mpi+mp)*(mpi+mp)))/(2*(mpi*mpi)))*18.0;           //WOLF's effective SIGMA
double RATIO_sig_bar = sigma_quasiel_PLUTO/sigma_quasiel_PLUTO1;  // RATIO OD ABOVE 2 effective sigma
//float sigma_quasiel_PLUTO= 2.*sigma_el*((s/(4.*mp*mp)) -1.);
//float dsigma=A*(sigma_quasiel/(Invmass*InvEn*InvEn))*(R2s2/R2s)*y*pt*pt;
float dsigma_PLUTO=A*(sigma_quasiel_PLUTO/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass);   //-----SCHAFER's FORMULA (dSigma/dM)
float dsigma_PLUTO_ph=A*(sigma_quasiel_PLUTO/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass)*(R2s2_M/R2s);  //SCHAFER's FORMULA (dSigma/dM)+PHASESPACE
float dsigma_PLUTO1=A*(sigma_quasiel_PLUTO/(Invmass*Invmass))*log((sqrts-(mpi+mp))/Invmass);  ////-----SCHAFER's FORMULA (dSigma/dM2)
float dsigma_PLUTO_Method2=WA2*(sigma_quasiel_PLUTO/(Invmass*pow(InvEnCM,3)))*(R2s2/R2s);  ////-----SCHAFER's FORMULA WORK_IN_PROGRESS (dSigma/d3p)  */
