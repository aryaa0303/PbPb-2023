import FWCore.ParameterSet.Config as cms

defaultCumu = cms.EDAnalyzer('MultiParticleCumulants', #Analyzer named: Correspond to the class name in 'plugin' folder
                             #Track collection
                            # tracks    = cms.InputTag('packedGenParticles'),
                             tracks    = cms.InputTag('packedPFCandidates'), 
                             tracksgen = cms.InputTag('packedGenParticles'), 
                             trackschi2 = cms.InputTag('packedPFCandidateTrackChi2'),
                             #qualityString = cms.string("highPurity"),
                             # tracks    = cms.InputTag('genParticles'),
                             #Vertex collection
                             #vertex    = cms.InputTag('offlinePrimaryVertices'),
                             vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             #                      qualityString = cms.string("highPurity"),
                             #Calorimeter tower collection
                            #caloTower = cms.InputTag('towerMaker'),
                             #Centrality
                             #entralitySrc    = cms.InputTag("hiCentrality"),
                             #entralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
                             centrality    = cms.InputTag("hiCentrality","","reRECO"),
                             centralitybin = cms.InputTag("centralityBin", "HFtowers"),
                             #Vertex selection
                             minvz         = cms.untracked.double(-15.0), 
                             maxvz         = cms.untracked.double(15.0),
                             maxrho        = cms.untracked.double(0.2),
                             isBVselByMult = cms.untracked.bool(True),
                             #Multiplicity selection
                             centmin       = cms.untracked.int32(0),  
                             centmax       = cms.untracked.int32(100),
                             noffmin       = cms.untracked.int32(10),
                             noffmax       = cms.untracked.int32(4000),
                             ptnoffmin     = cms.untracked.double(0.4),
                             ptnoffmax     = cms.untracked.double(10000.0),
                             dzdzerrornoff = cms.untracked.double(3.0),
                             d0d0errornoff = cms.untracked.double(3.0),
                             pterrorptnoff = cms.untracked.double(0.1),
                             #Track selection
                             etamin    = cms.untracked.double(-2.4),
                             etamax    = cms.untracked.double(2.4),
                             ptmin     = cms.untracked.double(0.5),
                             ptmax     = cms.untracked.double(3.0),
                             dzdzerror = cms.untracked.double(3.0),
                             d0d0error = cms.untracked.double(3.0),
                             pterrorpt = cms.untracked.double(0.1),
                             #Cumulant
                             cweight   = cms.untracked.bool(True),
                             branchSave = cms.untracked.int32(0),
                             #Acc X Eff
                             fname = cms.untracked.InputTag('General_TrackEff_3D_nhits0.root'),
                             effcentbin = cms.untracked.vint32(0,20,60,100,140,200),
                             fpt = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_highPt.root"),
                             fpix = cms.untracked.InputTag("2018PbPb_Efficiency_PixelTracks.root"),
                             fmb = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_MB.root"),
                             fplus = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_MB_ChargePlus.root"),
                             fminus = cms.untracked.InputTag("2018PbPb_Efficiency_GeneralTracks_MB_ChargeMinus.root")

)
