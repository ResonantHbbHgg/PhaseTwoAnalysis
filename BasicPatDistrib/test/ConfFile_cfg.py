
import FWCore.ParameterSet.Config as cms

process = cms.Process("MyAna")

PU = 0

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v1', '')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('MyAna')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('HH125GeV_200PU_3000.root') 
)



# 200 PU

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
 
  fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/02911C1E-9F72-E711-989D-A0369FC52064.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/06E5E1CE-3973-E711-9D15-001E67792882.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/0A1E9028-8272-E711-A01F-001E675A6AB3.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/0AE7E047-8072-E711-9AA5-008CFA110AD0.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/0E84FB2F-5072-E711-B2C9-549F35AD8BAF.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/10655306-2F73-E711-885F-7CD30AC030A2.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/120EDCB1-5F72-E711-AE82-0242AC130002.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/128FC9DD-0073-E711-A04C-002590D9D9E4.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/16D23AD6-5374-E711-AE07-B083FED0FFCF.root', 
###'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/PU200_91X_upgrade2023_realistic_v3-v3/110000/207B4B82-F972-E711-8F5E-60EB69BACA58.root'
)


)


if (PU == 0):

   process.source.fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/0A3BAD95-9E72-E711-BF10-002590D8C7BE.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/10F55CA1-4A73-E711-B109-001CC4A63C2A.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/1A973A4B-EF72-E711-9351-48FD8EE73AD1.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/2AD0A876-1173-E711-ABEC-3417EBE64567.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/3A7F1DBA-E072-E711-B7FC-002590E7D7E2.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/4098509E-4A73-E711-8B4A-141877642F9D.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/46083F80-4A73-E711-8D0A-549F3525A184.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/4654CE61-0B73-E711-BE91-90B11C2AA430.root',
###'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/529F8DE7-8D72-E711-92AF-0CC47A7C3424.root'
###'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring17MiniAOD/GluGluToHHTo2B2G_node_SM_14TeV-madgraph/MINIAODSIM/noPU_91X_upgrade2023_realistic_v3-v3/110000/5C1CED28-BB72-E711-A5DB-008CFA1112A0.root'

   )

   process.TFileService = cms.Service("TFileService",
        fileName = cms.string('HH125GeV_noPU_3000.root') 
   )





process.myana = cms.EDAnalyzer('BasicPatDistrib')
process.load("PhaseTwoAnalysis.BasicPatDistrib.CfiFile_cfi")
process.myana.useDeepCSV = True


process.p = cms.Path(process.myana)
