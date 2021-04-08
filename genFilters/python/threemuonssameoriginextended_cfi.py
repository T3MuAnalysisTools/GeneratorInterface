import FWCore.ParameterSet.Config as cms

threemufilter = cms.EDFilter("ThreeMuonsSameOriginExtended",
                                                NumRequired     = cms.int32(3),
                                                ParticleID      = cms.vint32(13,13,13),
                                                PtMin           = cms.vdouble(2.9,2.9,1.9),
                                                EtaMax          = cms.vdouble(2.45,2.45,2.45),
                                                Status          = cms.vint32(1,1,1),
                                                invMassMin      = cms.double(0),
                                                invMassMax      = cms.double(1000),
                                                maxDr           = cms.double(1.0)
                                            )




