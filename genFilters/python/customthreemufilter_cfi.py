import FWCore.ParameterSet.Config as cms

Custom3MuFilter = cms.EDFilter("CustomThreeMuFilter",
                               NumRequired     = cms.int32(3),
                               AcceptMore      = cms.bool(True),
                               ParticleID      = cms.vint32(13,13,13),
                               PtMin           = cms.vdouble(2.9,2.9,1.9),
                               EtaMax          = cms.vdouble(2.45,2.45,2.45),
                               Status          = cms.vint32(1,1,1),
                               invMassMin      = cms.double(1.5),
                               invMassMax      = cms.double(2.1)
                               )


