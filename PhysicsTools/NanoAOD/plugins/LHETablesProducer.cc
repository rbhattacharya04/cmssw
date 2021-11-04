#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/NanoAOD/interface/UniqueString.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "TLorentzVector.h"

#include <vector>
#include <iostream>
#include <boost/algorithm/string.hpp>

class LHETablesProducer : public edm::global::EDProducer<edm::BeginRunProducer> {
public:
  LHETablesProducer(edm::ParameterSet const& params)
      : lheLabel_(params.getParameter<std::vector<edm::InputTag>>("lheInfo")),
        lheTag_(edm::vector_transform(lheLabel_,
                                      [this](const edm::InputTag& tag) { return mayConsume<LHEEventProduct>(tag); })),
        lheRunTag_(edm::vector_transform(
            lheLabel_, [this](const edm::InputTag& tag) { return mayConsume<LHERunInfoProduct, edm::InRun>(tag); })),
        precision_(params.getParameter<int>("precision")),
        storeLHEParticles_(params.getParameter<bool>("storeLHEParticles")),
        storeAllLHEInfo_(params.getParameter<bool>("storeAllLHEInfo")) {
    produces<nanoaod::FlatTable>("LHE");
    if (storeLHEParticles_)
      produces<nanoaod::FlatTable>("LHEPart");
    if (storeAllLHEInfo_)
      produces<nanoaod::UniqueString,edm::InRun>("LHEInit");
  }

  ~LHETablesProducer() override {}
  
  void globalBeginRunProduce(edm::Run& iRun, edm::EventSetup const&) const override { 
    edm::Handle<LHERunInfoProduct> lheInfo;

    for (const auto& lheLabel : lheLabel_) {
      iRun.getByLabel(lheLabel, lheInfo);
      if (lheInfo.isValid()) {
        break;
      }
    }
    if (storeAllLHEInfo_ && lheInfo.isValid()) {
      // store LHE HEPRUP init
      auto hr = lheInfo->heprup();
      std::ostringstream strs;
      strs << hr.IDBMUP.first  << " ";
      strs << hr.IDBMUP.second << " ";
      strs << hr.EBMUP.first   << " ";
      strs << hr.EBMUP.second  << " ";
      strs << hr.PDFGUP.first  << " ";
      strs << hr.PDFGUP.second << " ";
      strs << hr.PDFSUP.first  << " ";
      strs << hr.PDFSUP.second << " ";
      strs << hr.IDWTUP        << " ";
      strs << hr.NPRUP         << " ";
      for (int i = 0; i < hr.NPRUP; ++i) {
        strs << "\n";
        strs << hr.XSECUP[i] << " ";
        strs << hr.XERRUP[i] << " ";
        strs << hr.XMAXUP[i] << " ";
        strs << hr.LPRUP[i]  << " ";
      }
    
      iRun.put(std::make_unique<nanoaod::UniqueString>(strs.str()), "LHEInit");
    }
  }

  void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {
    auto lheTab = std::make_unique<nanoaod::FlatTable>(1, "LHE", true);

    edm::Handle<LHEEventProduct> lheInfo;
    for (const auto& lheTag : lheTag_) {
      iEvent.getByToken(lheTag, lheInfo);
      if (lheInfo.isValid()) {
        break;
      }
    }
    if (lheInfo.isValid()) {
      auto lhePartTab = fillLHEObjectTable(*lheInfo, *lheTab);
      if (storeLHEParticles_)
        iEvent.put(std::move(lhePartTab), "LHEPart");
    } else {
      if (storeLHEParticles_) {  // need to store a dummy table anyway to make the framework happy
        auto lhePartTab = std::make_unique<nanoaod::FlatTable>(1, "LHEPart", true);
        iEvent.put(std::move(lhePartTab), "LHEPart");
      }
    }
    iEvent.put(std::move(lheTab), "LHE");
  }

  std::unique_ptr<nanoaod::FlatTable> fillLHEObjectTable(const LHEEventProduct& lheProd,
                                                         nanoaod::FlatTable& out) const {
    double lheHT = 0, lheHTIncoming = 0;
    unsigned int lheNj = 0, lheNb = 0, lheNc = 0, lheNuds = 0, lheNglu = 0;
    double lheVpt = 0;
    double alphaS = 0;
    double alphaQED=0;
    double scale  = 0;
    int    idproc = 0;
    int    rwl_type   = 0;
    int    rwl_index  = 0;
    int    rwl_seed   = 0;
    int    rwl_n1     = 0;
    int    rwl_n2     = 0;
    double rwl_weight = 0;

    const auto& hepeup = lheProd.hepeup();
    const auto& pup = hepeup.PUP;
    int lep = -1, lepBar = -1, nu = -1, nuBar = -1;
    std::vector<float> vals_pt;
    std::vector<float> vals_eta;
    std::vector<float> vals_phi;
    std::vector<float> vals_mass;
    std::vector<float> vals_pz;
    std::vector<float> vals_time;
    std::vector<int> vals_pid;
    std::vector<int> vals_status;
    std::vector<int> vals_spin;
    std::vector<int> vals_col1;
    std::vector<int> vals_col2;
    std::vector<int> vals_mother1;
    std::vector<int> vals_mother2;
    alphaS = hepeup.AQCDUP;
    alphaQED = hepeup.AQEDUP; 
    scale  = hepeup.SCALUP;
    idproc = hepeup.IDPRUP;
    if (storeAllLHEInfo_){
      for (unsigned int i = 0; i < lheProd.comments_size(); ++i) {
        std::string comment(lheProd.getComment(i));
        if (comment.find("#rwgt") != std::string::npos) {
          std::vector<std::string> vcomment;
          boost::algorithm::split(vcomment, comment, boost::is_any_of("\t "), boost::token_compress_on);
          if (vcomment.size() == 7) { // Powheg rwgt info
            rwl_type   = std::stoi(vcomment[1]);
            rwl_index  = std::stoi(vcomment[2]);
            rwl_weight = std::stod(vcomment[3]);
            rwl_seed   = std::stoi(vcomment[4]);
            rwl_n1     = std::stoi(vcomment[5]);
            rwl_n2     = std::stoi(vcomment[6]);
            break;
          }
        }
      }
    }
    for (unsigned int i = 0, n = pup.size(); i < n; ++i) {
      int status = hepeup.ISTUP[i];
      int idabs = std::abs(hepeup.IDUP[i]);
      if (status == 1 || status == -1 || storeAllLHEInfo_) {
        TLorentzVector p4(pup[i][0], pup[i][1], pup[i][2], pup[i][3]);  // x,y,z,t
        vals_pid.push_back(hepeup.IDUP[i]);
        vals_spin.push_back(hepeup.SPINUP[i]);
        vals_status.push_back(status);
        if (storeAllLHEInfo_){
            vals_col1.push_back(hepeup.ICOLUP[i].first);
            vals_col2.push_back(hepeup.ICOLUP[i].second);
            vals_mother1.push_back(hepeup.MOTHUP[i].first);
            vals_mother2.push_back(hepeup.MOTHUP[i].second);
            vals_time.push_back(hepeup.VTIMUP[i]);
        }
        if (status == -1) {
          vals_pt.push_back(0);
          vals_eta.push_back(0);
          vals_phi.push_back(0);
          vals_mass.push_back(0);
          vals_pz.push_back(p4.Pz());
        } else {
          vals_pt.push_back(p4.Pt());
          vals_eta.push_back(p4.Eta());
          vals_phi.push_back(p4.Phi());
          vals_mass.push_back(p4.M());
          vals_pz.push_back(0);
        }
      }
      if ((status == 1) && ((idabs == 21) || (idabs > 0 && idabs < 7))) {  //# gluons and quarks
        // object counters
        lheNj++;
        if (idabs == 5)
          lheNb++;
        else if (idabs == 4)
          lheNc++;
        else if (idabs <= 3)
          lheNuds++;
        else if (idabs == 21)
          lheNglu++;
        // HT
        double pt = std::hypot(pup[i][0], pup[i][1]);  // first entry is px, second py
        lheHT += pt;
        int mothIdx = std::max(
            hepeup.MOTHUP[i].first - 1,
            0);  //first and last mother as pair; first entry has index 1 in LHE; incoming particles return motherindex 0
        int mothIdxTwo = std::max(hepeup.MOTHUP[i].second - 1, 0);
        int mothStatus = hepeup.ISTUP[mothIdx];
        int mothStatusTwo = hepeup.ISTUP[mothIdxTwo];
        bool hasIncomingAsMother = mothStatus < 0 || mothStatusTwo < 0;
        if (hasIncomingAsMother)
          lheHTIncoming += pt;
      } else if (idabs == 12 || idabs == 14 || idabs == 16) {
        (hepeup.IDUP[i] > 0 ? nu : nuBar) = i;
      } else if (idabs == 11 || idabs == 13 || idabs == 15) {
        (hepeup.IDUP[i] > 0 ? lep : lepBar) = i;
      }
    }
    std::pair<int, int> v(0, 0);
    if (lep != -1 && lepBar != -1)
      v = std::make_pair(lep, lepBar);
    else if (lep != -1 && nuBar != -1)
      v = std::make_pair(lep, nuBar);
    else if (nu != -1 && lepBar != -1)
      v = std::make_pair(nu, lepBar);
    else if (nu != -1 && nuBar != -1)
      v = std::make_pair(nu, nuBar);
    if (v.first != -1 && v.second != -1) {
      lheVpt = std::hypot(pup[v.first][0] + pup[v.second][0], pup[v.first][1] + pup[v.second][1]);
    }

    out.addColumnValue<uint8_t>(
        "Njets", lheNj, "Number of jets (partons) at LHE step", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<uint8_t>("Nb", lheNb, "Number of b partons at LHE step", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<uint8_t>("Nc", lheNc, "Number of c partons at LHE step", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<uint8_t>(
        "Nuds", lheNuds, "Number of u,d,s partons at LHE step", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<uint8_t>(
        "Nglu", lheNglu, "Number of gluon partons at LHE step", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<float>("HT", lheHT, "HT, scalar sum of parton pTs at LHE step", nanoaod::FlatTable::FloatColumn);
    out.addColumnValue<float>("HTIncoming",
                              lheHTIncoming,
                              "HT, scalar sum of parton pTs at LHE step, restricted to partons",
                              nanoaod::FlatTable::FloatColumn);
    out.addColumnValue<float>("Vpt", lheVpt, "pT of the W or Z boson at LHE step", nanoaod::FlatTable::FloatColumn);
    out.addColumnValue<uint8_t>("NpNLO", lheProd.npNLO(), "number of partons at NLO", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<uint8_t>("NpLO", lheProd.npLO(), "number of partons at LO", nanoaod::FlatTable::UInt8Column);
    out.addColumnValue<float>("AlphaS", alphaS, "Per-event alphaS", nanoaod::FlatTable::FloatColumn);
    if (storeAllLHEInfo_){
        out.addColumnValue<float>("AlphaQED", alphaQED, "Per-event alphaQED", nanoaod::FlatTable::FloatColumn);
        out.addColumnValue<float>("Scale", scale, "Per-event scale", nanoaod::FlatTable::FloatColumn);
        out.addColumnValue<uint8_t>("ProcessID", idproc, "Process id (as in the card ordering)", nanoaod::FlatTable::UInt8Column);
        out.addColumnValue<int>("rwl_type", rwl_type,   "Per-event rwl_type for Powheg reweighting", nanoaod::FlatTable::IntColumn);
        out.addColumnValue<int>("rwl_index", rwl_index,  "Per-event rwl_index for Powheg reweighting", nanoaod::FlatTable::IntColumn);
        out.addColumnValue<float>("rwl_weight", rwl_weight, "Per-event rwl_weight for Powheg reweighting", nanoaod::FlatTable::FloatColumn);
        out.addColumnValue<int>("rwl_seed", rwl_seed,   "Per-event rwl_seed for Powheg reweighting", nanoaod::FlatTable::IntColumn);
        out.addColumnValue<int>("rwl_n1", rwl_n1,     "Per-event rwl_n1 for Powheg reweighting", nanoaod::FlatTable::IntColumn);
        out.addColumnValue<int>("rwl_n2", rwl_n2,     "Per-event rwl_n2 for Powheg reweighting", nanoaod::FlatTable::IntColumn);
    }

    auto outPart = std::make_unique<nanoaod::FlatTable>(vals_pt.size(), "LHEPart", false);
    outPart->addColumn<float>("pt", vals_pt, "Pt of LHE particles", nanoaod::FlatTable::FloatColumn, this->precision_);
    outPart->addColumn<float>(
        "eta", vals_eta, "Pseodorapidity of LHE particles", nanoaod::FlatTable::FloatColumn, this->precision_);
    outPart->addColumn<float>(
        "phi", vals_phi, "Phi of LHE particles", nanoaod::FlatTable::FloatColumn, this->precision_);
    outPart->addColumn<float>(
        "mass", vals_mass, "Mass of LHE particles", nanoaod::FlatTable::FloatColumn, this->precision_);
    outPart->addColumn<float>(
        "incomingpz", vals_pz, "Pz of incoming LHE particles", nanoaod::FlatTable::FloatColumn, this->precision_);
    outPart->addColumn<int>("pdgId", vals_pid, "PDG ID of LHE particles", nanoaod::FlatTable::IntColumn);
    outPart->addColumn<int>(
        "status", vals_status, "LHE particle status; -1:incoming, 1:outgoing", nanoaod::FlatTable::IntColumn);
    outPart->addColumn<int>("spin", vals_spin, "Spin of LHE particles", nanoaod::FlatTable::IntColumn);
    if (storeAllLHEInfo_){
        outPart->addColumn<int>("color1", vals_col1, "First color index of LHE particles", nanoaod::FlatTable::IntColumn);
        outPart->addColumn<int>("color2", vals_col2, "Second color index of LHE particles", nanoaod::FlatTable::IntColumn);
        outPart->addColumn<int>("mother1", vals_mother1, "First mother index of LHE particles", nanoaod::FlatTable::IntColumn);
        outPart->addColumn<int>("mother2", vals_mother2, "Second mother index of LHE particles", nanoaod::FlatTable::IntColumn);
        outPart->addColumn<float>("lifetime", vals_time, "Own lifetime of LHE particles", nanoaod::FlatTable::FloatColumn, this->precision_);
    }

    return outPart;
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::vector<edm::InputTag>>("lheInfo", std::vector<edm::InputTag>{{"externalLHEProducer"}, {"source"}})
        ->setComment("tag(s) for the LHE information (LHEEventProduct)");
    desc.add<int>("precision", -1)->setComment("precision on the 4-momenta of the LHE particles");
    desc.add<bool>("storeLHEParticles", false)
        ->setComment("Whether we want to store the 4-momenta of the status 1 particles at LHE level");
    desc.add<bool>("storeAllLHEInfo", false)->setComment("Whether to store the whole set of intermediate LHE particles, not only the status +/-1 ones");
    descriptions.add("lheInfoTable", desc);
  }

protected:
  const std::vector<edm::InputTag> lheLabel_;
  const std::vector<edm::EDGetTokenT<LHEEventProduct>> lheTag_;
  const std::vector<edm::EDGetTokenT<LHERunInfoProduct>> lheRunTag_;
  const unsigned int precision_;
  const bool storeLHEParticles_;
  const bool storeAllLHEInfo_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LHETablesProducer);
