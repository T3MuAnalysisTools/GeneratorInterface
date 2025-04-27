import FWCore.ParameterSet.Config as cms

Custom3PartFilter = cms.EDFilter("CustomThreeParticleFilter",
                               NumRequired     = cms.int32(2),
                               ParticleID      = cms.vint32(13,13),
                               PtMin           = cms.vdouble(1.0,1.0),
                               EtaMax          = cms.vdouble(2.9,2.9),
                               Status          = cms.vint32(1,1),
                               invMassMin      = cms.double(1.39),
                               invMassMax      = cms.double(2.11),
                               maxDr           = cms.double(1.2)
                               )


