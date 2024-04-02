// system include files                                                                                                                                                                                     
#include <memory>
#include <vector>

// CMSSW include files                                                                                                                                                                                      
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Analyzers/MultiParticleCumulants/data/EFF/trackingEfficiency2018PbPb.h"

// user include files                                                                                                                                                                                       
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include <TComplex.h>
#include <TObject.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
#include <TRandom3.h>

// class declaration                                                                                                                                                                                        
//                                                                                                                                                                                                          

// If the analyzer does not use TFileService, please remove                                                                                                                                                 
// the template argument to the base class so the class inherits                                                                                                                                            
// from  edm::one::EDAnalyzer<> and also remove the line from                                                                                                                                               
// constructor "usesResource("TFileService");"                                                                                                                                                              
// This will improve performance in multithreaded jobs.        

class TList;
class TF1;
class TH1;
class TH2;
class TH3F;
class TH1F;
class TH2F;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;


//class MultiParticleCumulants                                                                                                                                                                              
class MultiParticleCumulants : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  //member functions                                                                                                                                                                                          
 public:
  explicit MultiParticleCumulants(const edm::ParameterSet&);
  ~MultiParticleCumulants();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  int getEffNoffIndex();

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void FillQVector(TComplex _Qvector[20][20], double _Qcos[20][20], double _Qsin[20][20]);
  TComplex Two(int n1, int n2);
  //  TComplex TwoGap0p8(int n1, int n2);
  //TComplex TwoGap1(int n1, int n2);
  //TComplex TwoGap1p2(int n1, int n2);
  //TComplex TwoGap1p6(int n1, int n2);
  //TComplex TwoGap1p8(int n1, int n2);
  TComplex TwoGap2(int n1, int n2);
  TComplex Three(int n1, int n2, int n3);
  TComplex ThreeGap2_subM(int n1, int n2, int n3);
  TComplex ThreeGap2_subP(int n1, int n2, int n3);
  TComplex FourGap2(int n1, int n2, int n3, int n4);
  TComplex Four(int n1, int n2, int n3, int n4);
  TComplex Five(int n1, int n2, int n3, int n4, int n5);
  TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
  TComplex SixGap2(int n1, int n2, int n3, int n4, int n5, int n6) ;
  TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7) ; 
  TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
  TComplex Q(int n, int p);
  TComplex QGap0p8M(int n, int p);
  TComplex QGap0p8P(int n, int p);
  TComplex QGap1p0M(int n, int p);
  TComplex QGap1p0P(int n, int p);
  TComplex QGap1p2M(int n, int p);
  TComplex QGap1p2P(int n, int p);
  TComplex QGap1p6M(int n, int p);
  TComplex QGap1p6P(int n, int p);
  TComplex QGap1p8M(int n, int p);
  TComplex QGap1p8P(int n, int p);
  TComplex QGap2M(int n, int p);
  TComplex QGap2P(int n, int p);

  //  void ResetQ(const int nMaxHarm, const int nMaxPow);                                                                                                                                                   

  // ----------member data ---------------------------                                                                                                                                                      
  // ## tracks ##                                                                                                                                                                                           
  // used to select what tracks to read from configuration file                                                                                                                                             
  //edm::EDGetTokenT< edm::View < pat::PackedGenParticle > > trackTagsgen_; 
  //  edm::EDGetTokenT<reco::TrackCollection> trackTags_;
  //edm::EDGetTokenT<reco::GenParticleCollection> trackTagsgen_;                                                                                                                                         
  edm::EDGetTokenT< edm::View < pat::PackedCandidate > > trackTags_ ; 
  edm::EDGetTokenT< edm::View < pat::PackedGenParticle > > trackTagsgen_;
  edm::EDGetTokenT< edm::ValueMap < float > > chi2Map_;

  // ## vertex ##                                                                                                                                                                                           
  // used to select what vertex to read from configuration file                                                                                                                                             
  edm::EDGetTokenT<reco::VertexCollection> vtxTags_;
  // ## centrality ##                                                                                                                                                                                       
  // used to select what centrality collection to read from configuration file                                                                                                                              
  //edm::EDGetTokenT<reco::Centrality> centralityTags_;
  // used to access centrality bins                                                                                                                                                                        
  edm::EDGetTokenT<int> cent_bin_;


  //  edm::EDGetTokenT<int> centralityBinTags_;

  int cent_;

  // ## multiplicity selection (Noff)                                                                                                                                                                       
  int centmin_;
  int centmax_;
  int noffmin_;          //minimum multiplicity of an event to be considered                                                                                                                                
  int noffmax_;          //maximum multiplicity of an event to be considered                                                                                                                                
  double ptnoffmin_;     //minimum pt cut to compute Noff                                                                                                                                                   
  double ptnoffmax_;     //maximum pt cut to compute Noff                                                                                                                                                   
  double dzdzerrornoff_; //cut on dz/dzerror of the tracks to compute Noff                                                                                                                                  
  double d0d0errornoff_; //cut on d0/d0error of the tracks to compute Noff                                                                                                                                  
  double pterrorptnoff_; //cut on pterror/pt of the tracks to compute Noff                                                                                                                                  
  int noff_;             //ntrk offline value for a given event                                                                                                                                             

  // ## track selection ##                                                                                                                                                                                  
  double etamin_;    //min eta of the tracks                                                                                                                                                                
  double etamax_;    //max eta of the tracks                                                                                                                                                                
  double ptmin_;     //min pt of the tracks                                                                                                                                                                 
  double ptmax_;     //max pt of the tracks                                                                                                                                                                 
  double dzdzerror_; //cut on dz/dzerror of the tracks                                                                                                                                                      
  double d0d0error_; //cut on d0/d0error of the tracks                                                                                                                                                      
  double pterrorpt_; //cut on pterror/pt of the tracks                                                                                                                                                      
  int mult_;         //multiplicity (Nref) in a given event
  double mult_corr ;
  double pt_sum ; 

  // ## vertex selection ##                                                                                                                                                                                 
  double  minvz_;         //minimum z distance wrt (0,0,0) for the vertex                                                                                                                                   
  double  maxvz_;         //maximum z distance wrt (0,0,0) for the vertex                                                                                                                                   
  double  maxrho_;        //cut on rho distance for the vertex position                                                                                                                                     
  bool    isBVselByMult_; //sel best vertex based on vertex multiplicity (true) or sum p_T^2 (false)                                                                                                        
  int     nvtx_;          //number of reconstructed vertices in a given events                                                                                                                              
  double  xBestVtx_;      //x coordinate of the best vertex                                                                                                                                                 
  double  yBestVtx_;      //y coordinate of the best vertex                                                                                                                                                 
  double  rhoBestVtx_;    //r coordinate of the best vertex                                                                                                                                                 
  double  zBestVtx_;      //z coordinate of the best vertex                                                                                                                                                 
  double  xBestVtxError_; //x coordinate error of the best vertex                                                                                                                                           
  double  yBestVtxError_; //y coordinate error of the best vertex                                                                                                                                           
  double  zBestVtxError_; //z coordinate error of the best vertex                                                                                                                                           

  // ## harmonic and cumulants ##                                                                                                                                                                           
  bool cweight_; //use particle weight to correct from acc X eff                                                                                                                                            
  int branchSave_;  // whether force saving more branches                                                                                                                                                   

  // ## file acc & eff & fake ##                                                                                                                                                                            
  edm::InputTag fname_;         //file name that contains acc X eff corrections                                                                                                                             
  std::vector<int> effcentbin_ ; //Centrality binning of the correction                                                                                                                                     

  TFile* feff_;                 //TFile that contains 2D histos (pt, eta) with eff/(1-fak)                                                                                                                  
  //std::vector<TH2D*> heff_;     //TH2D histograms used for correction                   
  TH3D* hEff_3D;
  TH3D* hFak3D;

  //efficiency                                                                                                                                                                                         
  edm::InputTag fpt_;
  edm::InputTag fmb_;
  edm::InputTag fplus_;
  edm::InputTag fminus_;
  edm::InputTag fpix_;

  TrkEff2018PbPb* TrkEff;
  TrkEff2018PbPb* TrkEff1;
  TrkEff2018PbPb* TrkEff2;

  TH1F* hXBestVtx_;
  TH1F* hYBestVtx_;
  TH1F* hRhoBestVtx_;
  TH1F* hZBestVtx_;
  TH1F* hEtaTrk_;
  TH1F* hPtTrk_;
  TH1F* hEtaTrk_w;
  TH1F* hPtTrk_w;
  TH1F* hPhiTrk_;
  TH1F* hPhiTrk_w ; 
  TH1F* hEtaTrk_gen;
  TH1F* hPtTrk_gen;
  TH1F* hPhiTrk_gen; 
  TH1F* hEtaNoff_;
  TH1F* hPtNoff_;
  TH1F* hPhiNoff_;

  TH1F* V1_0_10;
  TH1F* V2_0_10;
  TH1F* S1_0_10;
  TH1F* n_in_Nch_0_10;

  TH1F* V1_0_5;
  TH1F* V2_0_5;
  TH1F* S1_0_5;
  TH1F* n_in_Nch_0_5;

  TH1F* V1_0_100;
  TH1F* V2_0_100;
  TH1F* S1_0_100;
  TH1F* n_in_Nch_0_100;
  
  // ## ttree ##                                                                                                                                                                                           

  //<<2>>                       
  TProfile* mean_pt ; 
 TProfile* fChcn22 ;
 TProfile* fChcn23 ;
 TProfile* fChcn24 ;
  TProfile* fChcn25 ;
  TProfile* fChcn26 ; 
  TProfile* wt2 ;
  
  TProfile* fChcn22_gap_0p8 ;
  TProfile* fChcn23_gap_0p8 ;
  TProfile* fChcn24_gap_0p8 ;
  TProfile* fChcn25_gap_0p8 ;
  TProfile* fChcn26_gap_0p8 ;
  TProfile* fChcn22_gap_1p0 ;
  TProfile* fChcn23_gap_1p0 ;
  TProfile* fChcn24_gap_1p0 ;
  TProfile* fChcn25_gap_1p0 ;
  TProfile* fChcn26_gap_1p0 ;
  TProfile* fChcn22_gap_1p2 ;
  TProfile* fChcn23_gap_1p2 ;
  TProfile* fChcn24_gap_1p2 ;
  TProfile* fChcn25_gap_1p2 ;
  TProfile* fChcn26_gap_1p2 ;
  TProfile* fChcn22_gap_1p6 ;
  TProfile* fChcn23_gap_1p6 ;
  TProfile* fChcn24_gap_1p6 ;
  TProfile* fChcn25_gap_1p6 ;
  TProfile* fChcn26_gap_1p6 ;
  TProfile* fChcn22_gap_1p8 ;
  TProfile* fChcn23_gap_1p8 ;
  TProfile* fChcn24_gap_1p8 ;
  TProfile* fChcn25_gap_1p8 ;
  TProfile* fChcn26_gap_1p8 ;
  TProfile* fChcn22_gap ;
  TProfile* fChcn23_gap ;
  TProfile* fChcn24_gap ;
  TProfile* fChcn25_gap ;
  TProfile* fChcn26_gap ;

  TProfile* wt2_gap_0p8 ;
  TProfile* wt2_gap_1p0 ;
  TProfile* wt2_gap_1p2 ;
  TProfile* wt2_gap_1p6 ;
  TProfile* wt2_gap_1p8 ;
  TProfile* wt2_gap ;   
  /*

 TProfile* fChcn22_gap ;
 TProfile* fChcn23_gap ;
 TProfile* fChcn24_gap ;
   TProfile* fChcn25_gap ;
  TProfile* fChcn26_gap ;
 TProfile* wt2_gap ;
  */
 //<<4>>                      

 TProfile* fChcn4_2 ;
 TProfile* fChcn4_3 ;
 TProfile* fChcn4_4 ;
 TProfile* wt4 ;

 TProfile* fChcn4_2_gap ;
 TProfile* fChcn4_3_gap ;
 TProfile* fChcn4_4_gap ;
 TProfile* wt4_gap ;

 TProfile* fChcnsc42 ;
 TProfile* fChcnsc32 ;
 TProfile* fChcnsc34 ;
 TProfile* fChcnsc25 ;
 TProfile* fChcnsc35 ;
 TProfile* fChcnsc45 ;
 TProfile* fChcnsc46 ;
 TProfile* fChcnsc26 ;

  //<<6>>                 
  
  TProfile* fChcn62 ;                                                                                                                                                                                                       
  TProfile*  fChcn63 ;                                                                                                                                                                                                       
  TProfile* fChcn64 ;                                                                                                                                                                                                       
  TProfile* fChcnmhc223 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc224 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc233 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc244 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc234 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc235 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc345 ;                                                                                                                                                                                                   
  TProfile* fChcnmhc246 ;                                                                                                                                                                                                   
  TProfile* wt6 ;                                                                                                                                                                                                           
                                                                                                                                                                                                                         
  TProfile* fChcn62_gap ;                                                                                                                                                                                                   
  TProfile* fChcn63_gap ;                                                                                                                                                                                                   
  TProfile* fChcn64_gap ;                                                                                                                                                                                                   
  TProfile* wt6_gap ;                                                                                                                                                                                                       
                                                                                                                                                                                                                         
  //<<8>>                                                                                                                                                                                                                
                                                                                                                                                                                                                         
  TProfile* fChcn28 ;                                                                                                                                                                                                       
  TProfile* fChcn38 ;                                                                                                                                                                                                       
  TProfile* fChcn48 ;                                                                                                                                                                                                       
  TProfile* fChcnmhc2223 ;                                                                                                                                                                                                  
  TProfile* fChcnmhc2333 ;                                                                                                                                                                                                  
  TProfile* fChcnmhc2224 ;                                                                                                                                                                                                  
  TProfile* fChcnmhc2444 ;                                                                                                                                                                                                  
  TProfile* fChcnmhc2233 ;                                                                                                                                                                                                  
  TProfile* fChcnmhc2244 ;                                                                                                                                                                                                  
  TProfile* wt8 ;                      
  
  //<<2>>	
  /*
  double fChcn22 ;
  double fChcn23 ;
  double fChcn24 ;
  double fChcn25 ;
  double fChcn26 ; 
  double wt2 ;


 double fChcn22_gap ;
 double fChcn23_gap ;
 double fChcn24_gap ;
 double fChcn25_gap ;
 double fChcn26_gap ;

 double wt2_gap ;

 //<<4>>                      

 double fChcn4_2 ;
 double fChcn4_3 ;
 double fChcn4_4 ;
 double wt4 ;

 double fChcn4_2_gap ;
 double fChcn4_3_gap ;
 double fChcn4_4_gap ;
 double wt4_gap ;

 double fChcnsc42 ;
 double fChcnsc32 ;
 double fChcnsc34 ;

 double fChcnsc42_gap ;
 double fChcnsc32_gap ;
 double fChcnsc34_gap ;


  //<<6>>                 
  
  double fChcn62 ;
  double fChcn63 ;
  double fChcn64 ; 
  double fChcnmhc223 ;
  double fChcnmhc224 ;
  double fChcnmhc233 ;
  double fChcnmhc244 ;
  double fChcnmhc234 ;
  double fChcnmhc235 ;
  double fChcnmhc345 ;
  double fChcnmhc246 ;
  double wt6 ;

  double fChcn62_gap ;
  double fChcn63_gap ;
  double fChcn64_gap ;
  double wt6_gap ;

  //<<8>>

  double fChcn28 ;
  double fChcn38 ;
  double fChcn48 ; 
  double fChcnmhc2223 ;
  double fChcnmhc2333 ;
  double fChcnmhc2224 ;
  double fChcnmhc2444 ;
  double fChcnmhc2233 ;
  double fChcnmhc2244 ;
  double wt8 ; 
  

   TProfile* meanpt_0_1 ;
  TProfile* meanpt_0_2 ;
  TProfile* meanpt_0_5 ;
  TProfile* meanpt_0_10 ;
  TProfile* meanpt_30_35 ;
  TProfile* meanpt_60_65 ;
  TProfile* meanpt_0_100 ;

  TProfile* meanpt_corr_0_1 ;
  TProfile* meanpt_corr_0_2 ;
  TProfile* meanpt_corr_0_5 ;
  TProfile* meanpt_corr_0_10 ;
  TProfile* meanpt_corr_30_35 ;
  TProfile* meanpt_corr_60_65 ;
  TProfile* meanpt_corr_0_100 ;

  TProfile* Nch_vs_cent_eta_m0p5_0p5;
  TProfile* meanpt_vs_cent_eta_m0p5_0p5;
  TProfile* Nch_vs_cent_eta_m2_2;
  TProfile* meanpt_vs_cent_eta_m2_2;
  
  TProfile* fChcn22_0_1 ;
  TProfile* fChcn22_0_2 ;
  TProfile* fChcn22_0_5 ;
  TProfile* fChcn22_0_10 ;
  TProfile* fChcn22_30_35;
  TProfile* fChcn22_60_65 ;
  TProfile* fChcn22_0_100 ;

  TProfile* fChcn4_2_0_1 ;
  TProfile* fChcn4_2_0_2 ;
  TProfile* fChcn4_2_0_5 ;
  TProfile* fChcn4_2_0_10 ;
  TProfile* fChcn4_2_30_35 ;
  TProfile* fChcn4_2_60_65 ;
  TProfile* fChcn4_2_0_100 ;

  TProfile* fChcn22_gap_0_1 ;
  TProfile* fChcn22_gap_0_2     ;
  TProfile* fChcn22_gap_0_5 ;
  TProfile* fChcn22_gap_0_10 ;
  TProfile* fChcn22_gap_30_35;
  TProfile* fChcn22_gap_60_65 ;
  TProfile* fChcn22_gap_0_100 ;
  */
  TTree* trEvent_;
  TTree* trEvent_gap ;
  edm::Service<TFileService> fs;
  
  TProfile* fmultvscent ;
  //TProfile* wt4_ ;
  //TProfile* wt6_ ;
  //TTree* trEvent_;
  //edm::Service<TFileService> fs;


  TComplex Qvector[20][20];

  TComplex Qvector0p8P[20][20] ;
  TComplex Qvector0p8M[20][20] ;
  TComplex Qvector1p0P[20][20] ;
  TComplex Qvector1p0M[20][20] ;
  TComplex Qvector1p2P[20][20] ;
  TComplex Qvector1p2M[20][20] ;
  TComplex Qvector1p6P[20][20] ;
  TComplex Qvector1p6M[20][20] ;
  TComplex Qvector1p8P[20][20] ;
  TComplex Qvector1p8M[20][20] ;
  TComplex Qvector2P[20][20] ;
  TComplex Qvector2M[20][20] ;
  double xbins[10000+10] = {}; //!                                                                                                                                                                         
  int nn = 0; //!                                                                                                                                                                                           

  // ClassDef(MultiParticleCumulants, 1) ;                                                                                                                                                                  

};

