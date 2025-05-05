import FWCore.ParameterSet.Config as cms
from pdfrecalculator_cfi import PDFRecalculator


# def addSVJPDFsAndScales(process,debug_flag=False, recalculatePDFs_flag=False, recalculateScales_flag=False, nEM_flag=0, nQCD_flag=0, pythiaSettings_flag='', pdfSetName_flag=''):
    
#     process.PDFRecalculator = PDFRecalculator.clone(
#             debug = cms.bool(debug_flag),
#             genEventInfoTag=cms.InputTag("generator"),
#             lheEventInfoTag=cms.InputTag("externalLHEProducer"),
#             recalculatePDFs = cms.bool(recalculatePDFs_flag),
#             recalculateScales = cms.bool(recalculateScales_flag),
#             pdfSetName = cms.string(pdfSetName_flag),
#         )
    
#     process.PDFRecalculator.nEM = nEM_flag
#     process.PDFRecalculator.recalculateScales = recalculateScales_flag
#     process.theoryWeightsMC = cms.Sequence(process.PDFRecalculator)

#     variables = {
#         "LHEPdfWeight" : ("LHEPdfWeight", "Float"),
#         "LHEScaleWeight" : ("LHEScaleWeight", "Float"),
#     }

#     for variable, (branch_name, type_) in variables.items():
#         process_name = variable + "Table"  # the module name MUST end by "Table"
#         setattr(process, process_name, cms.EDProducer(type_ + "ArrayTableProducer",
#             name = cms.string(branch_name),
#             doc = cms.string(""),
#             src = cms.InputTag("PDFRecalculator:" + variable),
#             )
#         )
#         process.theoryWeightsMC += getattr(process, process_name)

    
#     return process
def addSVJPDFsAndScales(process, debug_flag=False, recalculatePDFs_flag=False, recalculateScales_flag=False, nEM_flag=0, nQCD_flag=0, pythiaSettings_flag='', pdfSetName_flag=''):
    
    process.PDFRecalculator = PDFRecalculator.clone(
        debug = cms.bool(debug_flag),
        genEventInfoTag = cms.InputTag("generator"),
        lheEventInfoTag = cms.InputTag("externalLHEProducer"),
        recalculatePDFs = cms.bool(recalculatePDFs_flag),
        recalculateScales = cms.bool(recalculateScales_flag),
        pdfSetName = cms.string(pdfSetName_flag),
    )
    
    process.PDFRecalculator.nEM = nEM_flag
    process.PDFRecalculator.recalculateScales = recalculateScales_flag

    # Create a list to hold all modules for the sequence
    theory_weights_modules = [process.PDFRecalculator]

    variables = {
        "LHEPdfWeight": ("LHEPdfWeight", "Float"),
        "LHEScaleWeight": ("LHEScaleWeight", "Float"),
    }

    for variable, (branch_name, type_) in variables.items():
        process_name = variable + "Table"  # the module name MUST end by "Table"
        setattr(process, process_name, cms.EDProducer(type_ + "ArrayTableProducer",
            name = cms.string(branch_name),
            doc = cms.string(""),
            src = cms.InputTag("PDFRecalculator:" + variable),
        ))
        # Add the newly created module to the list
        theory_weights_modules.append(getattr(process, process_name))

    # Combine all modules into a single sequence using the '+' operator
    process.theoryWeightsMC = cms.Sequence(theory_weights_modules[0])
    for module in theory_weights_modules[1:]:
        process.theoryWeightsMC += module

    return process