// -*- C++ -*-
//
// Package:    Analyzers/Cumulants
// Class:      Cumulants
// 
/**\class Cumulants Cumulants.cc Analyzers/Cumulants/plugins/Cumulants.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maxime Guilbaud
//         Created:  Thu, 01 Jun 2017 16:56:11 GMT
//
//

// system include files

// CMSSW include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "Analyzers/MultiParticleCumulants/interface/MultiParticleCumulants.h"

//
// constructors and destructor
//
MultiParticleCumulants::MultiParticleCumulants(const edm::ParameterSet& iConfig) :

  //tracks & vertex
  // trackTags_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  trackTags_(consumes< edm::View< pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  trackTagsgen_(consumes< edm::View< pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("tracksgen"))),
  //  trackTagsgen_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  chi2Map_( consumes< edm::ValueMap< float > >( iConfig.getParameter< edm::InputTag >( "trackschi2" ) ) ),
  vtxTags_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  //centrality bin                                                                                                                                                                                  
  cent_bin_(consumes<int>(iConfig.getParameter<edm::InputTag>("centralitybin"))),
  //multiplicity selection
  centmin_(iConfig.getUntrackedParameter<int>("centmin")),
  centmax_(iConfig.getUntrackedParameter<int>("centmax")),
  noffmin_(iConfig.getUntrackedParameter<int>("noffmin")),
  noffmax_(iConfig.getUntrackedParameter<int>("noffmax")),
  ptnoffmin_(iConfig.getUntrackedParameter<double>("ptnoffmin")),
  ptnoffmax_(iConfig.getUntrackedParameter<double>("ptnoffmax")),
  dzdzerrornoff_(iConfig.getUntrackedParameter<double>("dzdzerrornoff")),
  d0d0errornoff_(iConfig.getUntrackedParameter<double>("d0d0errornoff")),
  pterrorptnoff_(iConfig.getUntrackedParameter<double>("pterrorptnoff")),
//track selection
  etamin_(iConfig.getUntrackedParameter<double>("etamin")),
  etamax_(iConfig.getUntrackedParameter<double>("etamax")),
  ptmin_(iConfig.getUntrackedParameter<double>("ptmin")),
  ptmax_(iConfig.getUntrackedParameter<double>("ptmax")),
  dzdzerror_(iConfig.getUntrackedParameter<double>("dzdzerror")),
  d0d0error_(iConfig.getUntrackedParameter<double>("d0d0error")),
  pterrorpt_(iConfig.getUntrackedParameter<double>("pterrorpt")),
//vertex selection
  minvz_(iConfig.getUntrackedParameter<double>("minvz")),
  maxvz_(iConfig.getUntrackedParameter<double>("maxvz")),
  maxrho_(iConfig.getUntrackedParameter<double>("maxrho")),
  isBVselByMult_(iConfig.getUntrackedParameter<bool>("isBVselByMult")),
  nvtx_(iConfig.getUntrackedParameter<int>("nvtx",-1)),
  xBestVtx_(iConfig.getUntrackedParameter<double>("xVtx",-99999)),
  yBestVtx_(iConfig.getUntrackedParameter<double>("yVtx",-99999)),
  rhoBestVtx_(iConfig.getUntrackedParameter<double>("rhoVtx",-99999)),
  zBestVtx_(iConfig.getUntrackedParameter<double>("zVtx",-99999)),
  xBestVtxError_(iConfig.getUntrackedParameter<double>("xVtxError",-99999)),
  yBestVtxError_(iConfig.getUntrackedParameter<double>("yVtxError",-99999)),
  zBestVtxError_(iConfig.getUntrackedParameter<double>("zVtxError",-99999)),
//harmonic order 
  cweight_(iConfig.getUntrackedParameter<bool>("cweight")),
  branchSave_(iConfig.getUntrackedParameter<int>("branchSave")),
//2-sub event relative eta difference
//  deltaeta_(iConfig.getUntrackedParameter<double>("deltaeta")),
//file acc & eff & fake
  fname_(iConfig.getUntrackedParameter<edm::InputTag>("fname")),
  effcentbin_(iConfig.getUntrackedParameter< std::vector<int> >("effcentbin")),

//EFF correction
  fpt_(iConfig.getUntrackedParameter<edm::InputTag>("fpt")),
  fmb_(iConfig.getUntrackedParameter<edm::InputTag>("fmb")),
  fplus_(iConfig.getUntrackedParameter<edm::InputTag>("fplus")),
  fminus_(iConfig.getUntrackedParameter<edm::InputTag>("fminus")),
  fpix_(iConfig.getUntrackedParameter<edm::InputTag>("fpix"))


{
  TString filename(fname_.label().c_str());
  feff_ = 0x0;
  if(cweight_ && !filename.IsNull())
    {
      edm::FileInPath fip(Form("Analyzers/MultiParticleCumulants/data/EFF/%s",filename.Data()));
      feff_ = new TFile(fip.fullPath().c_str(),"READ");

      //heff_.resize(feff_->GetNkeys());
      //if(heff_.size() != effcentbin_.size() - 1)
      //{
      //  edm::LogWarning ("Inconsitent binning") << " Inconsistent binning for the acc X eff correction..."
      //						  << " You might have wrong setting here";
      //	}

      //     for(unsigned short ik = 0; ik < heff_.size(); ++ik)
      //{
	  //heff_[ik] = (TH2D*) feff_->Get(feff_->GetListOfKeys()->At(ik)->GetName());
      //}

      hEff_3D = (TH3D*)feff_->Get("hEff_3D");
      hFak3D = (TH3D*)feff_->Get("hFak_3D");
      
      feff_->Close();
    }

  //**************************For efficiency correction ******************************************************
  /*                                            
  TString f_PT(fpt_.label().c_str());
  edm::FileInPath f1(Form("Analyzers/MultiParticleCumulants/data/EFF/effHpT/%s",f_PT.Data()));

  TString f_MB(fmb_.label().c_str());
  edm::FileInPath f2(Form("Analyzers/MultiParticleCumulants/data/EFF/effMB/%s",f_MB.Data()));

  TString f_Plus(fplus_.label().c_str());
  edm::FileInPath f3(Form("Analyzers/MultiParticleCumulants/data/EFF/plus/%s",f_Plus.Data()));

  TString f_Minus(fminus_.label().c_str());
  edm::FileInPath f4(Form("Analyzers/MultiParticleCumulants/data/EFF/minus/%s",f_Minus.Data()));
   
  TString f_Pix(fpix_.label().c_str());
  edm::FileInPath f5(Form("Analyzers/MultiParticleCumulants/data/EFF/effPix/%s",f_Pix.Data()));


  TrkEff = new TrkEff2018PbPb(  "general", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  TrkEff1 = new TrkEff2018PbPb(  "generalMB+", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  TrkEff2 = new TrkEff2018PbPb(  "generalMB-", false, f1.fullPath(), f2.fullPath(), f3.fullPath(), f4.fullPath(), f5.fullPath());
  */ 
  
  //Ouptut
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  // Histograms
}

void MultiParticleCumulants::beginJob()
{  
  //nn = 410 ;
  //for (int i = 0; i <= 410; i++) {
  // xbins[i] = (i + 0.5);
  // }
  // for (int i = 1; i <= 56; i++) {                                            
  // xbins[200+i] = (50*i + 200 + 0.5);                                         
  //}                                                                           

  nn = 7000 ;                                                                
  for(int i = 0 ; i <= 7000 ; i++) {                                           
    xbins[i] = i+0.5;                                                             
  }                     

  
  TFileDirectory fVtxHist  = fs->mkdir("Vertex");
  hXBestVtx_   = fVtxHist.make<TH1F>("hXvtx", "", 80, -0.2, 0.2);
  hYBestVtx_   = fVtxHist.make<TH1F>("hYvtx", "", 80, -0.2, 0.2);
  hRhoBestVtx_ = fVtxHist.make<TH1F>("hRvtx", "", 80, -0.2, 0.2);
  hZBestVtx_   = fVtxHist.make<TH1F>("hZvtx", "", 60, -30., 30.);
  //hcent_bin    = fVtxHist.make<TH1F>("hcent_bin", "", 200, 0.0, 200.0);
  TFileDirectory fTrkHist  = fs->mkdir("Tracks");
  //hEtaTrk_ = fTrkHist.make<TH1F>("hEtatrk", "", 30, -3.,   3.);
  hEtaTrk_ = fTrkHist.make<TH1F>("hEtatrk", "", 16, -3.2, 3.2);
  hPtTrk_  = fTrkHist.make<TH1F>("hPttrk",  "", 130,  0.,  13.);
  // hEtaTrk_w = fTrkHist.make<TH1F>("hEtatrkw", "", 30, -3.,   3.);
  hEtaTrk_w = fTrkHist.make<TH1F>("hEtatrkw", "", 16, -3.2, 3.2);
  hPtTrk_w  = fTrkHist.make<TH1F>("hPttrkw",  "", 130,  0.,  13.);
  hPhiTrk_ = fTrkHist.make<TH1F>("hPhitrk", "", 64, -3.2,  3.2);
  hPhiTrk_w = fTrkHist.make<TH1F>("hPhitrkw", "", 64, -3.2,  3.2);
  hEtaTrk_gen = fTrkHist.make<TH1F>("hEtatrk_gen", "", 16, -3.2, 3.2);
  hPtTrk_gen  = fTrkHist.make<TH1F>("hPttrk_gen",  "", 130,  0.,  13.);
  hPhiTrk_gen = fTrkHist.make<TH1F>("hPhitrk_gen", "", 64, -3.2,  3.2);
  hEtaNoff_ = fTrkHist.make<TH1F>("hEtaNoff", "", 30, -3.,   3.);
  hPtNoff_  = fTrkHist.make<TH1F>("hPtNoff",  "", 100,  0.,  200.);
  hPhiNoff_ = fTrkHist.make<TH1F>("hPhiNoff", "", 64, -3.2,  3.2);
  //TTree
  trEvent_ = fs->make<TTree>("trEvent", "trEvent");
  trEvent_->Branch("centrality", &cent_, "centrality/I");
  trEvent_->Branch("Noff",       &noff_, "Noff/I");
  trEvent_->Branch("Mult",       &mult_, "Mult/I");
  trEvent_->Branch("Corr_Mult",       &mult_corr, "Corr_Mult/D");
  /*
  trEvent_gap = fs->make<TTree>("trEvent_gap", "trEvent_gap") ;
  trEvent_gap->Branch("centrality", &cent_, "centrality/I");
  trEvent_gap->Branch("Noff",       &noff_, "Noff/I");
  trEvent_gap->Branch("Mult",       &mult_, "Mult/I");

  trEvent_gap->Branch("C22_2_gap", &fChcn22_gap, "C22_2_gap/D") ;
  trEvent_gap->Branch("C33_2_gap", &fChcn23_gap, "C33_2_gap/D") ;
  trEvent_gap->Branch("C44_2_gap", &fChcn24_gap, "C44_2_gap/D") ;
  trEvent_gap->Branch("C55_2_gap", &fChcn25_gap, "C55_2_gap/D") ;
  trEvent_gap->Branch("C66_2_gap", &fChcn26_gap, "C66_2_gap/D") ;
  trEvent_gap->Branch("W2_gap", &wt2_gap, "W2_gap/D") ;
 
  trEvent_gap->Branch("C22_4_gap", &fChcn4_2_gap, "C22_4_gap/D") ;
  trEvent_gap->Branch("C33_4_gap", &fChcn4_3_gap, "C33_4_gap/D") ;
  trEvent_gap->Branch("C44_4_gap", &fChcn4_4_gap, "C44_4_gap/D") ;
  trEvent_gap->Branch("W4_gap", &wt4_gap, "W4_gap/D") ;

  trEvent_gap->Branch("SC24_gap", &fChcnsc42_gap, "SC24_gap/D");
  trEvent_gap->Branch("SC23_gap", &fChcnsc32_gap, "SC23_gap/D");
  trEvent_gap->Branch("SC34_gap", &fChcnsc34_gap, "SC34_gap/D");

  trEvent_gap->Branch("C22_6_gap", &fChcn62_gap, "C22_6_gap/D") ;
  trEvent_gap->Branch("C33_6_gap", &fChcn63_gap, "C33_6_gap/D") ;
  trEvent_gap->Branch("C44_6_gap", &fChcn64_gap, "C44_6_gap/D") ;
  trEvent_gap->Branch("W6_gap", &wt6_gap, "W6_gap/D") ;


  
 
  //TFileDirectory fProfiles = fs->mkdir("fProfiles") ;                                               
  
//<<2>>                                                                                             
  trEvent_->Branch("C22_2", &fChcn22, "C22_2/D") ;
  trEvent_->Branch("C33_2", &fChcn23, "C33_2/D") ;
  trEvent_->Branch("C44_2", &fChcn24, "C44_2/D") ;
  trEvent_->Branch("C55_2", &fChcn25, "C55_2/D") ;
  trEvent_->Branch("C66_2", &fChcn26, "C66_2/D") ;
  trEvent_->Branch("W2", &wt2, "W2/D") ;
  */
   //<2>
  
  TFileDirectory fProfiles = fs->mkdir("fProfiles") ;    
  mean_pt = fProfiles.make<TProfile>("meanpt", "",  nn, xbins);
  /* fChcn22 = fProfiles.make<TProfile>("fchc22", "<cos 2(phi1 - phi2)>; Cent",  nn, xbins);
  fChcn23 = fProfiles.make<TProfile>("fchc23", "<cos 3(phi1 - phi2)> ; Cent",  nn, xbins);
  fChcn24 = fProfiles.make<TProfile>("fchc24", "<cos 4(phi1 - phi2)> ; Cent",  nn, xbins);
  fChcn25 = fProfiles.make<TProfile>("fchc25", "<cos 5(phi1 - phi2)> ; Cent",  nn, xbins);
  fChcn26 = fProfiles.make<TProfile>("fchc26", "<cos 6(phi1 - phi2)> ; Cent",  nn, xbins);
  wt2 = fProfiles.make<TProfile>("wt2", "w2 ; Cent", nn, xbins) ;
  fmultvscent = fProfiles.make<TProfile>("mvsc", "mvsc" , nn, xbins) ; 
  
   fChcn22_gap_0p8 = fProfiles.make<TProfile>("fchc22_gap_0p8", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn23_gap_0p8 = fProfiles.make<TProfile>("fchc23_gap_0p8", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn24_gap_0p8 = fProfiles.make<TProfile>("fchc24_gap_0p8", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn25_gap_0p8 = fProfiles.make<TProfile>("fchc25_gap_0p8", "<cos 5(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn26_gap_0p8 = fProfiles.make<TProfile>("fchc26_gap_0p8", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);

  fChcn22_gap_1p0 = fProfiles.make<TProfile>("fchc22_gap_1p0", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn23_gap_1p0 = fProfiles.make<TProfile>("fchc23_gap_1p0", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn24_gap_1p0 = fProfiles.make<TProfile>("fchc24_gap_1p0", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn25_gap_1p0 = fProfiles.make<TProfile>("fchc25_gap_1p0", "<cos 5(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn26_gap_1p0 = fProfiles.make<TProfile>("fchc26_gap_1p0", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);

  fChcn22_gap_1p2 = fProfiles.make<TProfile>("fchc22_gap_1p2", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn23_gap_1p2 = fProfiles.make<TProfile>("fchc23_gap_1p2", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn24_gap_1p2 = fProfiles.make<TProfile>("fchc24_gap_1p2", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn25_gap_1p2 = fProfiles.make<TProfile>("fchc25_gap_1p2", "<cos 5(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn26_gap_1p2 = fProfiles.make<TProfile>("fchc26_gap_1p2", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);

  fChcn22_gap_1p6 = fProfiles.make<TProfile>("fchc22_gap_1p6", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn23_gap_1p6 = fProfiles.make<TProfile>("fchc23_gap_1p6", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn24_gap_1p6 = fProfiles.make<TProfile>("fchc24_gap_1p6", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn25_gap_1p6 = fProfiles.make<TProfile>("fchc25_gap_1p6", "<cos 5(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn26_gap_1p6 = fProfiles.make<TProfile>("fchc26_gap_1p6", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);

  fChcn22_gap_1p8 = fProfiles.make<TProfile>("fchc22_gap_1p8", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn23_gap_1p8 = fProfiles.make<TProfile>("fchc23_gap_1p8", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn24_gap_1p8 = fProfiles.make<TProfile>("fchc24_gap_1p8", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn25_gap_1p8 = fProfiles.make<TProfile>("fchc25_gap_1p8", "<cos 5(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn26_gap_1p8 = fProfiles.make<TProfile>("fchc26_gap_1p8", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  
  fChcn22_gap = fProfiles.make<TProfile>("fchc22_gap_2p0", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn23_gap = fProfiles.make<TProfile>("fchc23_gap_2p0", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn24_gap = fProfiles.make<TProfile>("fchc24_gap_2p0", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn25_gap = fProfiles.make<TProfile>("fchc25_gap_2p0", "<cos 5(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
  fChcn26_gap = fProfiles.make<TProfile>("fchc26_gap_2p0", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);

  // wt2_gap_0p8 = fProfiles.make<TProfile>("wt2_gap_0p8", "w2 with eta gap ; Cent", nn, xbins) ;
  //wt2_gap_1p0 = fProfiles.make<TProfile>("wt2_gap_1p0", "w2 with eta gap ; Cent", nn, xbins) ;
  //wt2_gap_1p2 = fProfiles.make<TProfile>("wt2_gap_1p2", "w2 with eta gap ; Cent", nn, xbins) ;
  //wt2_gap_1p6 = fProfiles.make<TProfile>("wt2_gap_1p6", "w2 with eta gap ; Cent", nn, xbins) ;
  //wt2_gap_1p8 = fProfiles.make<TProfile>("wt2_gap_1p8", "w2 with eta gap ; Cent", nn, xbins) ;
 wt2_gap = fProfiles.make<TProfile>("wt2_gap_2p0", "w2 with eta gap ; Cent", nn, xbins) ;
  
 //fChcn24_gap = fProfiles.make<TProfile>("fchc24_gap", "<cos 4(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
 // fChcn25_gap = fProfiles.make<TProfile>("fchc25_gap", "<cos 5(phi1 - phi2)> with eta gap  ; Cent",  nn, xbins);
 // fChcn26_gap = fProfiles.make<TProfile>("fchc26_gap", "<cos 6(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
 //wt2_gap = fProfiles.make<TProfile>("wt2_gap", "w2 with eta gap ; Cent", nn, xbins) ;*/

 // fChcn22_gap_2p0 = fProfiles.make<TProfile>("fchc22_gap_2p0", "<cos 2(phi1 - phi2)> with eta gap; Cent",  nn, xbins);                                                                                    
 //  fChcn23_gap_2p0 = fProfiles.make<TProfile>("fchc23_gap_2p0", "<cos 3(phi1 - phi2)> with eta gap; Cent",  nn, xbins);
//<<4>>
/*
  trEvent_->Branch("C22_4", &fChcn4_2, "C22_4/D") ;                                                                                                                                                       
  trEvent_->Branch("C33_4", &fChcn4_3, "C33_4/D") ;                                                                                                                                                       
  trEvent_->Branch("C44_4", &fChcn4_4, "C44_4/D") ;                                                                                                                                                       
  trEvent_->Branch("SC42_4", &fChcnsc42, "SC42_4/D") ;                                                                                                                                                    
  trEvent_->Branch("SC32_4", &fChcnsc32, "SC32_4/D") ;                                                                                                                                                    
  trEvent_->Branch("SC34_4", &fChcnsc34, "SC34/D") ; 
  trEvent_->Branch("W4", &wt4, "W4/D") ;                                                                                                                                                                   

  //<<6>>

  trEvent_->Branch("C22_6", &fChcn62, "C22_6/D") ;
  trEvent_->Branch("C33_6", &fChcn63, "C33_6/D") ;
  trEvent_->Branch("C44_6", &fChcn64, "C44_6/D") ;
  trEvent_->Branch("MHC_223", &fChcnmhc223, "MHC_223/D") ;                                                                                                                                                
  trEvent_->Branch("MHC_224", &fChcnmhc224, "MHC_224/D") ;
  trEvent_->Branch("MHC_233", &fChcnmhc233, "MHC_233/D") ;    
  trEvent_->Branch("MHC_244", &fChcnmhc244, "MHC_244/D") ;  
  trEvent_->Branch("MHC_234", &fChcnmhc234, "MHC_234/D") ;
  trEvent_->Branch("MHC_235", &fChcnmhc235, "MHC_235/D") ;
  trEvent_->Branch("MHC_345", &fChcnmhc345, "MHC_345/D") ;
  trEvent_->Branch("MHC_246", &fChcnmhc246, "MHC_246/D") ;
  trEvent_->Branch("W6", &wt6, "W6/D") ;                                                                                                                                                                   

  //<<8>>

  trEvent_->Branch("C22_8", &fChcn28, "C22_8/D") ;
  trEvent_->Branch("C33_8", &fChcn38, "C33_8/D") ;
  trEvent_->Branch("C44_8", &fChcn48, "C44_8/D") ;
  trEvent_->Branch("MHC_2223", &fChcnmhc2223, "MHC_2223/D") ;
  trEvent_->Branch("MHC_2333", &fChcnmhc2333, "MHC_2333/D") ;
  trEvent_->Branch("MHC_2224", &fChcnmhc2224, "MHC_2224/D") ;
  trEvent_->Branch("MHC_2444", &fChcnmhc2444, "MHC_2444/D") ;
  trEvent_->Branch("MHC_2233", &fChcnmhc2233, "MHC_2233/D") ;
  trEvent_->Branch("MHC_2244", &fChcnmhc2244, "MHC_2244/D") ;
  trEvent_->Branch("W8", &wt8, "W8/D") ; 

  
  //<4>                                                                                                                                                                                                     
  
  fChcn4_2 = fProfiles.make<TProfile>("fchc4_2", "<cos 2(phi1 + phi2 - phi3 - phi4)> ; Cent",  nn, xbins);
  fChcn4_3 = fProfiles.make<TProfile>("fchc4_3", "<cos 3(phi1 + phi2 - phi3 - phi4)> ; Cent",  nn, xbins);
  fChcn4_4 = fProfiles.make<TProfile>("fchc4_4", "<cos 4(phi1 + phi2 - phi3 - phi4)> ; Cent",  nn, xbins);
  wt4 = fProfiles.make<TProfile>("wt4", "w4 ; Centrality(%)", nn, xbins) ;

  fChcn4_2_gap = fProfiles.make<TProfile>("fchc4_2_gap", "<cos 2(phi1 + phi2 - phi3 - phi4)> with eta gap ; Cent",  nn, xbins);
  fChcn4_3_gap = fProfiles.make<TProfile>("fchc4_3_gap", "<cos 3(phi1 + phi2 - phi3 - phi4)> with eta gap; Cent",  nn, xbins);
  fChcn4_4_gap = fProfiles.make<TProfile>("fchc4_4_gap", "<cos 4(phi1 + phi2 - phi3 - phi4)> with eta gap; Cent",  nn, xbins);
  wt4_gap = fProfiles.make<TProfile>("wt4_gap", "w4 with eta gap ; Centrality(%)", nn, xbins) ;

  //sc                                                                                                                                                                                                      

  fChcnsc42 = fProfiles.make<TProfile>("fchcsc42", "<cos(2*phi1 + 4*phi2 - 2*phi3 - 4*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc32 = fProfiles.make<TProfile>("fchcsc32", "<cos(2*phi1 + 3*phi2 - 2*phi3 - 3*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc34 = fProfiles.make<TProfile>("fchcsc34", "<cos(3*phi1 + 4*phi2 - 3*phi3 - 4*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc25 = fProfiles.make<TProfile>("fchcsc25", "<cos(2*phi1 + 5*phi2 - 2*phi3 - 5*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc35 = fProfiles.make<TProfile>("fchcsc35", "<cos(3*phi1 + 5*phi2 - 3*phi3 - 5*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc45 = fProfiles.make<TProfile>("fchcsc45", "<cos(4*phi1 + 5*phi2 - 4*phi3 - 5*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc46 = fProfiles.make<TProfile>("fchcsc46", "<cos(4*phi1 + 6*phi2 - 4*phi3 - 6*phi4)> ; Cent(%)",  nn, xbins);
  fChcnsc26 = fProfiles.make<TProfile>("fchcsc26", "<cos(2*phi1 + 6*phi2 - 2*phi3 - 6*phi4)> ; Cent(%)",  nn, xbins) ;

   //<6>                                                                                                                                                                                                     
  //mhc(2,23)                                                                                                                                                                                              
  fChcnmhc223 = fProfiles.make<TProfile>("fchcmhc223", "<cos(2*phi1 + 2*phi2 + 3*phi3 - 2*phi4 - 2*phi5 - 3*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc224 = fProfiles.make<TProfile>("fchcmhc224", "<cos(2*phi1 + 2*phi2 + 4*phi3 - 2*phi4 - 2*phi5 - 4*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc233 = fProfiles.make<TProfile>("fchcmhc233", "<cos(2*phi1 + 3*phi2 + 3*phi3 - 2*phi4 - 3*phi5 - 3*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc244 = fProfiles.make<TProfile>("fchcmhc244", "<cos(2*phi1 + 4*phi2 + 4*phi3 - 2*phi4 - 4*phi5 - 4*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc234 = fProfiles.make<TProfile>("fchcmhc234", "<cos(2*phi1 + 3*phi2 + 4*phi3 - 2*phi4 - 3*phi5 - 4*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc235 = fProfiles.make<TProfile>("fchcmhc235", "<cos(2*phi1 + 3*phi2 + 5*phi3 - 2*phi4 - 3*phi5 - 5*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc345 = fProfiles.make<TProfile>("fchcmhc345", "<cos(3*phi1 + 4*phi2 + 5*phi3 - 3*phi4 - 4*phi5 - 5*phi6)> ; Cent(%)", nn, xbins);
  fChcnmhc246 = fProfiles.make<TProfile>("fchcmhc246", "<cos(2*phi1 + 4*phi2 + 6*phi3 - 2*phi4 - 4*phi5 - 6*phi6)> ; Cent(%)", nn, xbins);
  wt6 = fProfiles.make<TProfile>("wt6", "w6 ; Centrality(%)", nn, xbins) ;  

  fChcn62 = fProfiles.make<TProfile>("fchcn62", "<cos(2*phi1 + 2*phi2 + 2*phi3 - 2*phi4 - 2*phi5 - 2*phi6> ; Cent", nn, xbins);
  fChcn63 = fProfiles.make<TProfile>("fchcn63", "<cos(3*phi1 + 3*phi2 + 3*phi3 - 3*phi4 - 3*phi5 - 3*phi6> ; Cent", nn, xbins);
  fChcn64 = fProfiles.make<TProfile>("fchcn64", "<cos(4*phi1 + 4*phi2 + 4*phi3 - 4*phi4 - 4*phi5 - 4*phi6> ; Cent", nn, xbins);

  wt6_gap = fProfiles.make<TProfile>("wt6_gap", "w6_gap ; Centrality(%)", nn, xbins) ;
  fChcn62_gap = fProfiles.make<TProfile>("fchcn62_gap", "<cos(2*phi1 + 2*phi2 + 2*phi3 - 2*phi4 - 2*phi5 - 2*phi6> ; Cent", nn, xbins);
  fChcn63_gap = fProfiles.make<TProfile>("fchcn63_gap", "<cos(3*phi1 + 3*phi2 + 3*phi3 - 3*phi4 - 3*phi5 - 3*phi6> ; Cent", nn, xbins);
  fChcn64_gap = fProfiles.make<TProfile>("fchcn64_gap", "<cos(4*phi1 + 4*phi2 + 4*phi3 - 4*phi4 - 4*phi5 - 4*phi6> ; Cent", nn, xbins);

  //<8>
  fChcn28 = fProfiles.make<TProfile>("fchcn28", " ; cent", nn, xbins) ; 
  fChcn38 = fProfiles.make<TProfile>("fchcn38", " ; cent", nn, xbins) ;
  fChcn48 = fProfiles.make<TProfile>("fchcn48", " ; cent", nn, xbins) ;	
  wt8 = fProfiles.make<TProfile>("wt8", "w8 ; Centrality(%)", nn, xbins) ;

  fChcnmhc2223 = fProfiles.make<TProfile>("fchcnmhc2223", " ; cent", nn, xbins) ;
  fChcnmhc2233 = fProfiles.make<TProfile>("fchcnmhc2233", " ; cent", nn, xbins) ;
  fChcnmhc2333 = fProfiles.make<TProfile>("fchcnmhc2333", " ; cent", nn, xbins) ;
  fChcnmhc2224 = fProfiles.make<TProfile>("fchcnmhc2224", " ; cent", nn, xbins) ;
  fChcnmhc2244 = fProfiles.make<TProfile>("fchcnmhc2244", " ; cent", nn, xbins) ;
  fChcnmhc2444 = fProfiles.make<TProfile>("fchcnmhc2444", " ; cent", nn, xbins) ;
*/
 
}
   

MultiParticleCumulants::~MultiParticleCumulants()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MultiParticleCumulants::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
  using namespace std;

  //track collection 
  //edm::Handle< reco::TrackCollection > tracks;
  //iEvent.getByToken(trackTags_, tracks);
  auto trks = iEvent.getHandle( trackTags_ );
  auto trksgen = iEvent.getHandle( trackTagsgen_ ) ; 

  //edm::Handle< reco::GenParticleCollection > tracks;
  //iEvent.getByToken(trackTagsgen_, tracks);

  double NtrksAfterGap0p8M = 0 ;
  double NtrksAfterGap0p8P = 0 ;
  double NtrksAfterGap1p0M = 0 ;
  double NtrksAfterGap1p0P = 0 ;
  double NtrksAfterGap1p2M = 0 ;
  double NtrksAfterGap1p2P = 0 ;
  double NtrksAfterGap1p6M = 0 ;
  double NtrksAfterGap1p6P = 0 ;
  double NtrksAfterGap1p8M = 0 ;
  double NtrksAfterGap1p8P = 0 ;
  double NtrksAfterGap2p0M = 0 ;
  double NtrksAfterGap2p0P = 0 ;

  //standard correlation
  /*double Qcos[20][20] ;                                                                            
  double Qsin[20][20] ;                                                                                                                                                                                  
  //eta gap, from left                                                                                                                                                                                    
  /* double QcosGap0p8M[20][20] ;
  double QsinGap0p8M[20][20] ;
  double QcosGap1p0M[20][20] ;
  double QsinGap1p0M[20][20] ;
  double QcosGap1p2M[20][20] ;
  double QsinGap1p2M[20][20] ;
  double QcosGap1p6M[20][20] ;
  double QsinGap1p6M[20][20] ;
  double QcosGap1p8M[20][20] ;
  double QsinGap1p8M[20][20] ;
  
  double QcosGap2M[20][20] ;
  double QsinGap2M[20][20] ;

  //eta gap, from right                                                                                                                                                                                  
  /* double QcosGap0p8P[20][20] ;
  double QsinGap0p8P[20][20] ;
  double QcosGap1p0P[20][20] ;
  double QsinGap1p0P[20][20] ;
  double QcosGap1p2P[20][20] ;
  double QsinGap1p2P[20][20] ;
  double QcosGap1p6P[20][20] ;
  double QsinGap1p6P[20][20] ;
  double QcosGap1p8P[20][20] ;
  double QsinGap1p8P[20][20] ;
  
  double QcosGap2P[20][20] ;
  double QsinGap2P[20][20] ;
  
  for(int i = 0 ; i < 20 ; i++)
    { for(int j = 0 ; j < 20 ; j++)
	{
	  Qcos[i][j] = 0 ;                                                                        
	  Qsin[i][j] = 0 ;
	  /*  QcosGap0p8M[i][j] = 0 ;
          QsinGap0p8M[i][j] = 0 ;
          QcosGap0p8P[i][j] = 0 ;
          QsinGap0p8P[i][j] = 0 ;
	  QcosGap1p0M[i][j] = 0 ;
          QsinGap1p0M[i][j] = 0 ;
          QcosGap1p0P[i][j] = 0 ;
          QsinGap1p0P[i][j] = 0 ;
	  QcosGap1p2M[i][j] = 0 ;
          QsinGap1p2M[i][j] = 0 ;
          QcosGap1p2P[i][j] = 0 ;
          QsinGap1p2P[i][j] = 0 ;
	  QcosGap1p6M[i][j] = 0 ;
          QsinGap1p6M[i][j] = 0 ;
          QcosGap1p6P[i][j] = 0 ;
          QsinGap1p6P[i][j] = 0 ;
	  QcosGap1p8M[i][j] = 0 ;
          QsinGap1p8M[i][j] = 0 ;
          QcosGap1p8P[i][j] = 0 ;
          QsinGap1p8P[i][j] = 0 ;
	  
	  QcosGap2M[i][j] = 0 ;
	  QsinGap2M[i][j] = 0 ;
	  QcosGap2P[i][j] = 0 ;
	  QsinGap2P[i][j] = 0 ;
	  
	}
    }
  */
  
  double V1_term = 0;
  double V2_term = 0;
  double S1_term = 0;

  double Q1 = 0;
  double Q2 = 0;
  double Q3 = 0;
  
  //access tracks chi2/ndf
  auto chi2Map = iEvent.getHandle( chi2Map_ );

  //vtx collection      
  auto pvs = iEvent.getHandle( vtxTags_ );
  //best vertex                                                                                                                                                                                 
  //nvtx_  = 0;     //N valid vertex in collection                            
  //edm::Handle< reco::VertexCollection > vertices;
  //iEvent.getByToken(vtxTags_, vertices);
  // reco::VertexCollection verticesColl = *vertices;
  //VertexCollection vsorted = *vertices;
  // const reco::VertexCollection *recoV = verticesColl.product();
  
  //if( !vertices->size() )
  //{
  //  edm::LogWarning ("Missing Collection") <<"Invalid or empty vertex collection!";
  //  return;
  // }
 
               
  xBestVtx_      = -999.;
  yBestVtx_      = -999.;
  rhoBestVtx_    = -999.;
  zBestVtx_      = -999.; //Best vtx coordinates                                                     
  xBestVtxError_ = -999.;
  yBestVtxError_ = -999.;
  zBestVtxError_ = -999.; //Best vtx error      
         
  double bestvzError;
  math::XYZPoint bestvtx;
  math::Error<3>::type vtx_cov;
   if ( !pvs->empty() ) {
     const reco::Vertex& vtx = (*pvs)[0];
  //  int primaryvtx = 0 ; 
 
     bestvzError = vtx.zError();
     bestvtx = vtx.position();
     vtx_cov = vtx.covariance();
     }else {
     return;
   }

  xBestVtx_ = bestvtx.x();
  yBestVtx_ = bestvtx.y();
  zBestVtx_ = bestvtx.z();
   
 //  hZBestVtx_ -> Fill(zBestVtx_);
  
  //if ( zBestVtx_ < -15.0 || zBestVtx_ >= 15.0 ) return;
  if ( std::abs( bestvtx.z() ) > 15. ) return; //wide
  hZBestVtx_ -> Fill(zBestVtx_);
  
  // ----------------- centrality selection -------------------------------                                                                                                                                 //access centrality bins 
  auto cbin = iEvent.getHandle( cent_bin_ );
  double centBin = ( float ) (*cbin);
  
  if (centBin < 0 || centBin >= 10) return;
  cent_ = centBin ; 
  //hZBestVtx_ -> Fill(zBestVtx_);
  //hcent_bin -> Fill(centBin);

  
  //Reset QVectors to start fresh
  //qN_.reset();
  mult_ = 0; // Event multiplicity
  mult_corr = 0 ; 
  //double mult_small_ = 0;                                                                               
  std::vector<double> val(2,0.);

   pt_sum = 0 ;
   // double pt_sum_small = 0;
   //double pt_sum_corr = 0 ;
   //double wt_sum = 0 ;
   
   // auto cbin = iEvent.getHandle( cent_bin_ );                                                                                                                                                                
   // double centBin = ( float ) (*cbin);                                                                                                                                                                       
  
   //if (centBin < 0 || centBin >= 200) return;                                                                                                                                                                
   // cent_ = centBin ;      
   
  //noff_ = 0 ; 
  //track loop
   int trkIndx = -1;

   for(auto const& iter_tk : *trksgen)
     {
       trkIndx++;

       //if ( !trk.hasTrackDetails() ) continue;
       //auto iter_tk = trk.pseudoTrack();

       double pt = iter_tk.pt();
       double eta = iter_tk.eta();
       int charge = iter_tk.charge();
       double phi = iter_tk.phi();
       
       if( charge == 0 ) continue;
       if( pt < 0.0001) continue ;

       if(eta < -2.4 || eta > 2.4) continue ;
       if(pt < 0.3 || pt > 3.0)    continue ;
       if(iter_tk.status() != 1) continue;

       hEtaTrk_gen->Fill(eta);
       hPtTrk_gen ->Fill(pt);
       hPhiTrk_gen->Fill(phi);
     }

   for (auto const& trk : *trks)
 // for( reco::TrackCollection::const_iterator iter_tk = tracks->begin();
   //    iter_tk != tracks->end(); ++iter_tk )  
   // for( reco::GenParticleCollection::const_iterator itTrk = tracks->begin();
   //   itTrk != tracks->end();
   //   ++itTrk )
     {
       trkIndx++;

       if ( !trk.hasTrackDetails() ) continue;
       auto iter_tk = trk.pseudoTrack();
      
       double dzvtx = iter_tk.dz( bestvtx );
       double dxyvtx = iter_tk.dxy( bestvtx );
       double dzerror = std::hypot( iter_tk.dzError(), bestvzError );
       double dxyerror = iter_tk.dxyError( bestvtx, vtx_cov );
       double pterror = iter_tk.ptError();

       //double dzvtx = iter_tk->dz( bestvtx );
    //double dxyvtx = iter_tk->dxy( bestvtx );
    //double dzerror = std::hypot( iter_tk->dzError(), bestvzError );
    //double dxyerror = iter_tk->dxyError( bestvtx, vtx_cov );
       // double pterror = iter_tk->ptError();

      // Get eta, pt, and charge of the track                                                                                                                                                              
       double pt = iter_tk.pt();
       double eta = iter_tk.eta();
       int charge = iter_tk.charge();
       double phi = iter_tk.phi();

       //double pt = itTrk->pt();
       //double eta = itTrk->eta();
       //int charge = itTrk->charge();
       //double phi = itTrk->phi();

       auto hit_pattern = iter_tk.hitPattern();

//HI specific cuts                                                                                                                                                                                   
       double chi2ndof = ( double ) ( *chi2Map )[ trks->ptrAt( trkIndx ) ];
       double dcaxy = (dxyvtx / dxyerror);
       double dcaz = (dzvtx / dzerror);
       double ptreso = (fabs(pterror) / pt);
       int nHits = iter_tk.numberOfValidHits();
       double chi2n = ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() );

//selected tracks                                                                                                                                                                                    
    if( charge == 0 ) continue;
    if( pt < 0.0001) continue ;
    //if(itTrk->status() != 1) continue ; 
    if( fabs(pterror) / pt > 0.1 ) continue; //changed, 30.07.2023 
    if( fabs(dzvtx / dzerror) >= 3.0 ) continue; 
    if( fabs(dxyvtx / dxyerror) >= 3.0  ) continue; 
    if( ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() ) >= 0.18 ) continue;
    //if ( iter_tk.numberOfValidHits() < 11 ) continue;

    //tight and loose
    //if( fabs(pterror) / pt > 0.15 ) continue; //changed, 30.07.2023               
    //if( fabs(dzvtx / dzerror) >= 5.0 ) continue;                                     
    //if( fabs(dxyvtx / dxyerror) >= 5.0  ) continue;                                  
    //if( ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() ) >= 0.18 ) continue;
    //if ( iter_tk.numberOfValidHits() < 11 ) continue;                                
    if( iter_tk.quality(reco::TrackBase::qualityByName("highPurity")) != 1 ) continue;
	
    // Track selection for analysis                                                                                                                                                                    
    if(eta < -2.4 || eta > 2.4) continue ;
    if(pt < 0.3 || pt > 3.0)    continue ;
    
    pt_sum += pt ;

    // Compute weights
    //double weight = -999.0 ;
    double weight = 1.0 ; 

    double eff = hEff_3D->GetBinContent(hEff_3D->GetXaxis()->FindBin(eta), hEff_3D->GetYaxis()->FindBin(pt), hEff_3D->GetZaxis()->FindBin(centBin));
    double fake = hFak3D->GetBinContent(hEff_3D->GetXaxis()->FindBin(eta), hEff_3D->GetYaxis()->FindBin(pt), hEff_3D->GetZaxis()->FindBin(centBin));

    //std::cout << pt << ", " << eff << ", " << fake << std::endl ; 
    weight = (1 - fake)/eff;
    //std::cout << weight << std::endl ; 
    // if (weight == 0) {
    //std::cout << weight << std::endl ;
    // }
    
    //++mult_;
    mult_ += 1 ; 
    mult_corr += weight ; 
    //pt_sum += weight * pt ;
    //wt_sum += weight ;
  
    //   Q1 += pt;
    //Q2 += pow(pt,2);
    //Q3 += pow(pt,3);

    //if(eta > -0.5 && eta < 0.5)
    // {
    //	pt_sum_small += pt;
    //	++mult_small_;
    // }
    
    // Fill trk histograms
    // if (!std::isfinite(eta) || !std::isfinite(pt) || !std::isfinite(phi) || !std::isfinite(weight)) continue ; 
    hEtaTrk_w->Fill(eta,weight);
    hPtTrk_w ->Fill(pt,weight);
    hPhiTrk_w->Fill(phi,weight);
    hEtaTrk_->Fill(eta);
    hPtTrk_ ->Fill(pt);
    hPhiTrk_->Fill(phi);
     
    //..calculate Q-vectors
    //..no eta gap
    /*
    for(int iharm=0; iharm<16; iharm++) {
      for(int ipow=0; ipow<9; ipow++) {
    	Qcos[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
    	Qsin[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
      }
    }
       /*
    if(eta < -0.4) {
      NtrksAfterGap0p8M++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap0p8M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap0p8M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
    if(eta > 0.4) {
      NtrksAfterGap0p8P++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap0p8P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap0p8P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
	}
      }
    }

    if(eta < -0.5) {
      NtrksAfterGap1p0M++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p0M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p0M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
    if(eta > 0.5) {
      NtrksAfterGap1p0P++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p0P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p0P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
	}
      }
    }

    if(eta < -0.6) {
      NtrksAfterGap1p2M++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p2M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p2M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
    if(eta > 0.6) {
      NtrksAfterGap1p2P++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p2P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p2P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }

    if(eta < -0.8) {
      NtrksAfterGap1p6M++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p6M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p6M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
    if(eta > 0.8) {
      NtrksAfterGap1p6P++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p6P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p6P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }

    if(eta < -0.9) {
      NtrksAfterGap1p8M++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p8M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p8M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
    if(eta > 0.9) {
      NtrksAfterGap1p8P++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap1p8P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap1p8P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
       
    if(eta < -1) {
      NtrksAfterGap2p0M++;
      for(int iharm=0; iharm<16; iharm++) {
	for(int ipow=0; ipow<6; ipow++) {
	  QcosGap2M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
	  QsinGap2M[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
	}
      }
    }
    if(eta > 1) {
      NtrksAfterGap2p0P++;
      for(int iharm=0; iharm<16; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
	  QcosGap2P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*phi);
          QsinGap2P[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*phi);
        }
      }
    }
    */
    
     }
    
  //end of track loop

  // if(mult_corr == 0) return ;  
  // mean_pt->Fill(mult_corr, pt_sum/mult_corr) ; 
   //FillQVector(Qvector, Qcos, Qsin) ;
   //FillQVector(Qvector2M, QcosGap2M, QsinGap2M);                                                                                                                                                       
   //FillQVector(Qvector2P, QcosGap2P, QsinGap2P); 
   //  fmultvscent->Fill(cent_/2, mult_corr) ; 
   /*
   double Dn2 = 0 ;
   Dn2 =  Two(0, 0).Re();

   //standard two particle correlation                                                                                                                                                                      

     if(mult_ > 1 && Dn2 != 0) {

     //v2{2}                                                                                                                                                                                               
   
     TComplex v22 =  Two(2, -2);
     double v22Re = v22.Re()/Dn2;                                                                                                                                                                                          
     fChcn22->Fill(cent_/2, v22Re, Dn2);
     fChcn22->Sumw2() ;
     wt2->Fill(cent_/2, Dn2, Dn2);
  
     //v3{2}                                                                                                                                                                                               
     TComplex v32 =  Two(3, -3);
     double v32Re = v32.Re()/Dn2;
     fChcn23->Fill(cent_/2, v32Re, Dn2);

     //v4{2}                                                                                                                                                                                                
     TComplex v42 =  Two(4, -4);
     double v42Re = v42.Re()/Dn2;
     fChcn24->Fill(cent_/2, v42Re, Dn2);

     TComplex v52 =  Two(5, -5);
     double v52Re = v52.Re()/Dn2;                                                                                                                                                                        
     fChcn25->Fill(cent_/2, v52Re, Dn2);

     TComplex v62 =  Two(6, -6);
     double v62Re = v62.Re()/Dn2;                                                                                                                                                                        
     fChcn26->Fill(cent_/2, v62Re, Dn2); 

   }

   double Dn4 = 0 ;
   Dn4 = Four(0,0,0,0).Re() ;
   
   //standard four particle correlation                                                                                                                                                                     
   
   if(mult_ > 3 && Dn4 != 0) {
     //vn{4}                                                                                                                                                                                                
     
     TComplex v24 = Four(2, 2, -2, -2) ;
     double v24Re = v24.Re()/Dn4 ;
     fChcn4_2->Fill(cent_/2, v24Re, Dn4) ; 
     wt4->Fill(cent_/2, Dn4, Dn4);

     TComplex v34 = Four(3, 3, -3, -3) ;
     double v34Re = v34.Re()/Dn4 ;
     fChcn4_3->Fill(cent_/2, v34Re, Dn4);

     TComplex v44 = Four(4, 4, -4, -4) ;
     double v44Re = v44.Re()/Dn4 ;
     fChcn4_4->Fill(cent_/2, v44Re, Dn4);
      
     TComplex sc42 = Four(2, 4, -2, -4) ;
     double sc42Re = sc42.Re()/Dn4;
     fChcnsc42->Fill(cent_/2, sc42Re, Dn4) ;
      
     TComplex sc32 = Four(2, 3,-2, -3) ;
     double sc32Re = sc32.Re()/Dn4 ;
     fChcnsc32->Fill(cent_/2, sc32Re, Dn4) ;
    
     TComplex sc34 = Four(3, 4,-3, -4) ;
     double sc34Re = sc34.Re()/Dn4 ;
     fChcnsc34->Fill(cent_/2, sc34Re, Dn4) ;

     TComplex sc25 = Four(2, 5, -2, -5) ;
     double sc25Re = sc25.Re()/Dn4 ;
     fChcnsc25->Fill(cent_/2, sc25Re, Dn4) ;

     TComplex sc35 = Four(3, 5, -3, -5) ;
     double sc35Re = sc35.Re()/Dn4 ;
     fChcnsc35->Fill(cent_/2, sc35Re, Dn4) ;

     TComplex sc45 = Four(4, 5, -4, -5) ;
     double sc45Re = sc45.Re()/Dn4 ;
     fChcnsc45->Fill(cent_/2, sc45Re, Dn4) ;

     TComplex sc46 = Four(4, 6, -4, -6) ;
     double sc46Re = sc46.Re()/Dn4 ;
     fChcnsc46->Fill(cent_/2, sc46Re, Dn4) ;

     TComplex sc26 = Four(2, 6, -2, -6) ;
     double sc26Re = sc26.Re()/Dn4 ;
     fChcnsc26->Fill(cent_/2, sc26Re, Dn4) ;

   }

   double Dn6 = 0 ;
   Dn6 = Six(0,0,0,0,0,0).Re() ;
    
    //mixed harmonics                                                                                                                                                                                        
    
   if(mult_ > 5 && Dn6 != 0) {
      
     //mhc(2,2,3)                                                                                                                                                                                           
     TComplex mhc223 = Six(2, 2, 3, -2, -2, -3);
     double mhc223Re = mhc223.Re()/Dn6 ;
     fChcnmhc223->Fill(cent_/2, mhc223Re, Dn6) ;                                                                                                                                                                          
     wt6->Fill(cent_/2, Dn6, Dn6) ;                                                                                                                                                                                          
      
      //mhc(2,2,4)                                                                                                                                                                                           
     TComplex mhc224 = Six(2, 2, 4, -2, -2, -4);
     double mhc224Re = mhc224.Re()/Dn6 ;    
     fChcnmhc224->Fill(cent_/2, mhc224Re, Dn6);
     
     
     //mhc(2,3,3)                                                                                                                                                                                           
     TComplex mhc233 = Six(2, 3, 3, -2, -3, -3);
     double mhc233Re = mhc233.Re()/Dn6 ;
     //double mhc233Re = mhc233.Re() ;
     //fChcnmhc233 = mhc233Re ;    
     fChcnmhc233->Fill(cent_/2, mhc233Re, Dn6);
     
     //mhc(2,4,4)                                                                                                                                                                                           
     TComplex mhc244 = Six(2, 4, 4, -2, -4, -4);
     double mhc244Re = mhc244.Re()/Dn6 ;
     //double mhc244Re = mhc244.Re() ;
     //fChcnmhc244 = mhc244Re ;    
     fChcnmhc244->Fill(cent_/2, mhc244Re, Dn6);

       //mhc(2,3,4)                                                                                                                                                                                           
     TComplex mhc234 = Six(2, 3, 4, -2, -3, -4);
     double mhc234Re = mhc234.Re()/Dn6 ;
     //double mhc234Re = mhc234.Re() ; 
     //fChcnmhc234 = mhc234Re ;    
     fChcnmhc234->Fill(cent_/2, mhc234Re, Dn6);

     TComplex mhc235 = Six(2, 3, 5, -2, -3, -5);
     double mhc235Re = mhc235.Re()/Dn6 ;                                                                                                                                                               
     //double mhc235Re = mhc235.Re() ;
     fChcnmhc235->Fill(cent_/2, mhc235Re, Dn6) ;

     TComplex mhc345 = Six(3, 4, 5, -3, -4, -5);
     double mhc345Re = mhc345.Re()/Dn6 ;                                                                                                                                                               
     //double mhc345Re = mhc345.Re() ;
     //fChcnmhc345 = mhc345Re ;
     fChcnmhc345->Fill(cent_/2, mhc345Re, Dn6) ;
 
     TComplex mhc246 = Six(2, 4, 6, -2, -4, -6);
     double mhc246Re = mhc246.Re()/Dn6 ;                                                                                                                                                               
     //double mhc246Re = mhc246.Re() ;
     fChcnmhc246->Fill(cent_/2, mhc246Re, Dn6) ;

     TComplex v26 = Six(2, 2, 2, -2, -2, -2);
     double v26Re = v26.Re()/Dn6 ;
     fChcn62->Fill(cent_/2, v26Re, Dn6) ;

     TComplex v36 = Six(3, 3, 3, -3, -3, -3);
     double v36Re = v36.Re()/Dn6 ;
     fChcn63->Fill(cent_/2, v36Re, Dn6) ;

     TComplex v46 = Six(4, 4, 4, -4, -4, -4);
     double v46Re = v46.Re()/Dn6 ;
     fChcn64->Fill(cent_/2, v46Re, Dn6) ;


      
   }
                                                                                                                                                                                           
    double Dn8 = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
    
    if(mult_ > 7 && Dn8 != 0)
      {
	  TComplex v28 = Eight(2, 2, 2, 2, -2, -2, -2, -2);
	  double v28Re = v28.Re()/Dn8;
	  fChcn28->Fill(cent_/2, v28Re, Dn8) ;
	  wt8->Fill(cent_/2, Dn8, Dn8) ;

	  TComplex v38 = Eight(3, 3, 3, 3, -3, -3, -3, -3);
          double v38Re = v38.Re()/Dn8;
	  fChcn38->Fill(cent_/2, v38Re, Dn8) ;

	  TComplex v48 = Eight(4, 4, 4, 4, -4, -4, -4, -4);
          double v48Re = v48.Re()/Dn8;
	  fChcn48->Fill(cent_/2, v48Re, Dn8) ; ;

	  TComplex mhc2223 = Eight(2, 2, 2, 3, -2, -2, -2, -3);
	  double mhc2223Re = mhc2223.Re()/Dn8 ;
	  fChcnmhc2223->Fill(cent_/2, mhc2223Re, Dn8) ;

	  TComplex mhc2333 = Eight(2, 3, 3, 3, -2, -3, -3, -3);
          double mhc2333Re = mhc2333.Re()/Dn8 ;
          fChcnmhc2333->Fill(cent_/2, mhc2333Re, Dn8) ;

	  TComplex mhc2224 = Eight(2, 2, 2, 4, -2, -2, -2, -4);
          double mhc2224Re = mhc2224.Re()/Dn8 ;
	  fChcnmhc2224->Fill(cent_/2, mhc2224Re, Dn8) ;
	  
	  TComplex mhc2444 = Eight(2, 4, 4, 4, -2, -4, -4, -4);
	  double mhc2444Re = mhc2444.Re()/Dn8 ;
	  fChcnmhc2444->Fill(cent_/2, mhc2444Re, Dn8) ;
	  
	  TComplex mhc2233 = Eight(2, 2, 3, 3, -2, -2, -3, -3);
	  double mhc2233Re = mhc2233.Re()/Dn8 ;
	  fChcnmhc2233->Fill(cent_/2, mhc2233Re, Dn8) ;

	  TComplex mhc2244 = Eight(2, 2, 4, 4, -2, -2, -4, -4);
          double mhc2244Re = mhc2244.Re()/Dn8 ;
	  fChcnmhc2244->Fill(cent_/2, mhc2244Re, Dn8) ;

      }
  
   
    double Dn2Gap2 = 0 ;
   Dn2Gap2 = TwoGap2(0,0).Re() ;

   if(NtrksAfterGap2p0M > 0 && NtrksAfterGap2p0P > 0 && Dn2Gap2 != 0)
     {
       TComplex v22Gap2 = TwoGap2(2,-2) ;
       //double v22ReGap2 = v22Gap2.Re() ;                                                                                                       \
                                                                                                                                                   
       //fChcn22_gap = v22ReGap2 ;                                                                                                               \
                                                                                                                                                   
       //wt2_gap = Dn2Gap2 ;                                                                                                                     \
                                                                                                                                                   
       double v22ReGap2 = v22Gap2.Re()/Dn2Gap2 ;
       //fChcn22_mult -> Fill(mult_, v22ReGap2, Dn2Gap2) ;                                                                                        \
                                                                                                                                                   
       fChcn22_gap->Fill(cent_/2, v22ReGap2, Dn2Gap2) ;
       wt2_gap->Fill(cent_/2, Dn2Gap2, Dn2Gap2);

       TComplex v32Gap2 = TwoGap2(3,-3) ;
       //double v32ReGap2 = v32Gap2.Re() ;                                                                                                       \
                                                                                                                                                   
       //fChcn23_gap = v32ReGap2 ;                                                                                                                
       double v32ReGap2 = v32Gap2.Re()/Dn2Gap2 ;
       //fChcn23_mult -> Fill(mult_, v32ReGap2, Dn2Gap2) ;                                                                                        \
                                                                                                                                                   
       fChcn23_gap->Fill(cent_/2, v32ReGap2, Dn2Gap2) ;

       TComplex v42Gap2 = TwoGap2(4,-4) ;
       double v42ReGap2= v42Gap2.Re()/Dn2Gap2 ;                                                                                                  \

       fChcn24_gap->Fill(cent_/2, v42ReGap2, Dn2Gap2) ;                                                                                          \


       TComplex v52Gap2 = TwoGap2(5,-5) ;
       double v52ReGap2= v52Gap2.Re()/Dn2Gap2 ;
       fChcn25_gap->Fill(cent_/2, v52ReGap2, Dn2Gap2) ;

       TComplex v62Gap2 = TwoGap2(6,-6) ;
       double v62ReGap2= v62Gap2.Re()/Dn2Gap2 ;
       fChcn26_gap->Fill(cent_/2, v62ReGap2, Dn2Gap2) ;

     }
   */
    cent_ = centBin ;
    trEvent_->Fill() ;
   //trEvent_gap->Fill() ;
   
} //end of event loop

  // ------------ method called once each job just after ending the event loop  ------------
void 
MultiParticleCumulants::endJob() 
{
}

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MultiParticleCumulants::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int
MultiParticleCumulants::getEffNoffIndex() 
{
   for( int idx = 0; idx < (int) effcentbin_.size() - 1; ++idx )
     {
       if( cent_ >= effcentbin_[idx] && cent_ < effcentbin_[idx+1] ) return idx;
     }
   return -1;  
}

TComplex MultiParticleCumulants::Q(int n, int p)
{

  if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);
  
}

TComplex MultiParticleCumulants::QGap2M(int n, int p)
{

  if(n>=0) return Qvector2M[n][p];
  else return TComplex::Conjugate(Qvector2M[-n][p]);

}

TComplex MultiParticleCumulants::QGap2P(int n, int p)
{

  if(n>=0) return Qvector2P[n][p];
  else return TComplex::Conjugate(Qvector2P[-n][p]);

}

TComplex MultiParticleCumulants::Two(int n1, int n2)
{
  
  TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
  
}

TComplex MultiParticleCumulants::TwoGap2(int n1, int n2)
{

  TComplex formula = QGap2M(n1,1)*QGap2P(n2,1);
  return formula;

}

TComplex MultiParticleCumulants::Three(int n1, int n2, int n3)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
		- Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
	return formula;

}

TComplex MultiParticleCumulants::Four(int n1, int n2, int n3, int n4)
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
    + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
  return formula;
}

TComplex MultiParticleCumulants::Five(int n1, int n2, int n3, int n4, int n5)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
	return formula;

}

TComplex MultiParticleCumulants::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{
  
  
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
    - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
    - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
    + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
    + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
    - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
    - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
    - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
    - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
    + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
    + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
    - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
    - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
    - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
    - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
    + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
    - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
    - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
    - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
    - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
    + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
    - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
    - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
    - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
    + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
    - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
    + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
    - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
    - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
    - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
    + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
    + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
    - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
    - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
    - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
    - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
    + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
    + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
    - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
    - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
    - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
    - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
    + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
    - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
    - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
    - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
    + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
    + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
    + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
    + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
    + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
    + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
    + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
    - 120.*Q(n1+n2+n3+n4+n5+n6,6);
  return formula;
  
}

TComplex MultiParticleCumulants::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6};

	for(int k=7; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[6] = {0,1,2,3,4,5};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==6: there is just one combination, we can add it manually
		if(k==6){
		  Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
						Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])])*
						Q(Narray[int(array[5])]+n7, 7-k);
				}
			}while(std::next_permutation(array, array+6));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
		}// k==0

		else{
			std::cout<<"invalid range of k"<<std::endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;
}

TComplex MultiParticleCumulants::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6, n7};

	for(int k=8; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[7] = {0,1,2,3,4,5,6};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==7: there is just one combination, we can add it manually
		if(k==7){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
				Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
		}// k==7

		else if(k==6){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
						Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
						Q(Narray[int(array[6])]+n8, 8-k);
				}
			}while(std::next_permutation(array, array+7));
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])], Narray[int(array[4])])*
							Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%120 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
		}// k==0

		else{
			std::cout<<"invalid range of k"<<std::endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}

void MultiParticleCumulants::FillQVector(TComplex _Qvector[20][20], double _Qcos[20][20], double _Qsin[20][20]) {
  for(int iharm=0; iharm<20; iharm++)
    {
      for(int ipow=0; ipow<20; ipow++)
	{
	  _Qvector[iharm][ipow] = TComplex(_Qcos[iharm][ipow], _Qsin[iharm][ipow]);
	}
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(MultiParticleCumulants);
