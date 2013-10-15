#include "JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "JetMETObjects/interface/SimpleJetCorrector.h"
#include "JetMETObjects/interface/JetCorrectorParameters.h"
#include <vector>
 
JetCorrectorParameters corr;
JetCorrectorParameters::Definitions def;
JetCorrectorParameters::Record record;
std::vector<JetCorrectorParameters> corrv;
std::vector<JetCorrectorParameters::Record> recordv;
JetCorrectorParametersCollection coll;
JetCorrectorParametersCollection::pair_type pair_type;
JetCorrectorParametersCollection::collection_type colltype;
std::vector<JetCorrectorParametersCollection> collv;
