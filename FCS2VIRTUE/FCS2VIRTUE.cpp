// FCS2VIRTUE: converts a SimpleTree ROOT file's FCS calorimeter hits (and
// optionally mcparticle tracks and reco/truth jets from a corresponding
// jetTree file) into the JSON event format read by VIRTUE for event
// display/animation. C++ port of FCS2VIRTUE.py (same geometry constants,
// same per-hit block/hit classification and energy-color scale), reading
// trees directly via ROOT instead of uproot so it can run as a compiled
// standalone tool rather than needing a Python+uproot environment.
//
// Usage:
//   ./FCS2VIRTUE <input_simpletree.root> <n_events> <output.json> [options]
//
// Options:
//   --tracks           include mcparticle tracks (SimpleTree's mcpart_*
//                      branches)
//   --truth-jets       include truth jets (green), read from --jets-file
//   --reco-jets        include reco jets (orange), read from --jets-file
//   --jets-file <file> the jetTree file to read jets from; required if
//                      either --truth-jets or --reco-jets is given.
//                      Truth and reco jets always come from the same
//                      file (both live in the same jetTree), so this is
//                      a single shared path, not one per jet type.
//
// --jets-file must be the 1:1 JetMatcher/RecoJets output for THIS
// SimpleTree file (same event ordering) -- see data_to_jet's/
// sim_to_jet's suffix-matching convention. Omit --truth-jets/--reco-jets
// entirely if you don't want jets; no jetTree file is needed in that case.
//
// Compile (needs the fastjet module loaded, since ../shared/JetParameters.h
// pulls in FastJet headers even though this file doesn't cluster anything
// itself -- it only reads that header's shared R constant):
//   module load fastjet-3.3.4
//   g++ -m64 FCS2VIRTUE.cpp ../shared/VectorDict.cxx -o FCS2VIRTUE `root-config --cflags --libs` `fastjet-config --cxxflags --libs` -std=c++11
//
// ../shared/VectorDict.cxx/.h provide the compiled vector<float>/
// vector<int> CollectionProxy dictionary jetTree's branches need (same
// files jet_scale/JetEnergyScaleFineGrid.cpp links against to read the
// same jetTree schema); without it, TTree::SetBranchAddress on a
// vector<float> branch fails at runtime with "do not have a compiled
// CollectionProxy". Referenced from shared/ rather than duplicated
// here since this is a standalone local tool, not a SUMS pipeline that
// needs to be self-contained for SandBox packaging. Only needed for
// --truth-jets/--reco-jets; the SimpleTree hit/track branches are plain
// arrays and don't require it, but it's harmless to always link.
//
// ../shared/JetParameters.h provides R (the anti-kt clustering radius in
// eta-phi space, shared with JetMatcher.cpp/RecoJets.cpp) so this file
// never hardcodes its own copy that could drift out of sync with
// whatever radius actually produced the jetTree being rendered.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstring>

#include "TFile.h"
#include "TTree.h"

#include "../shared/JetParameters.h" // shared R (jet clustering radius, eta-phi)

// ============================================================
// Geometry constants (identical to FCS2VIRTUE.py)
// ============================================================

static const double UNIT_SCALE = 10.0; // ROOT units are cm -> output mm

static const double ECAL_SIZE[3] = {55.2, 55.2, 330.0};
static const double HCAL_SIZE[3] = {100.0, 100.0, 840.0};
static const double EPD_HIT_SIZE = 50.0;

static const double C_LIGHT = 299.792458; // mm/ns
static const double DETECTOR_ANGLE = 1.73003; // degrees

static const int MAX = 10000; // max FCS hits/mcparticles per event, same buffer convention as JetMatcher.cpp

// ------------------------------------------------------------
// Tracker/jet-length geometry, from FCS_Model.json (its length_unit is
// meters; everything below is pre-converted to mm to match this file's
// own length_unit).
//   - tracker_boundary: per VIRTUE's own docs, this is a cylindrical
//     boundary [radius, left_boundary, right_boundary] -- tracks
//     terminate once they cross it. Sourced from FCS_Model's "STAR Main
//     Detector" toroid (index 5): position [0,0,0], radii up to 2.20 m
//     (-> radius), length [4.50, 4.50] m -- each of those two values is
//     the TOTAL length on that side, not a half-length, so the actual
//     boundary offset from position is length/2 = 2.25 m (-> left/right
//     boundary at -2.25 m / +2.25 m).
//   - FCS_FRONT_Z_MM: FCS_Model's ECAL Left/Right blocks (index 0/1) sit
//     at position z = 7.398065 m with a 0.330 m depth along z; that
//     position is the block's CENTER (matching shift_block_position's
//     front-face -> center convention below), so the front face is at
//     center_z - depth/2 = 7.398065 - 0.165 = 7.233065 m. ECAL is the
//     closer of the two FCS calorimeters (HCAL sits further out, at
//     8.476607 m), so this is "the front FCS surface" a jet gets drawn
//     out to.
// ------------------------------------------------------------

static const double TRACKER_RADIUS_MM = 2200.0;
static const double TRACKER_Z_MIN_MM = -2250.0;
static const double TRACKER_Z_MAX_MM = 2250.0;

static const double FCS_FRONT_Z_MM = 7233.065;

// STAR's standard solenoid field strength, 0.5 T -- the simulation chain
// (runSimBfc.C's input_chain) uses the "FieldOn" BFC option, which is
// STAR's normal full-field configuration, not a reduced or zero field.
static const double B_FIELD_TESLA = 0.5;

// VIRTUE's track-propagation step size (time between successive points
// along a rendered helix).
static const double SEGMENT_NS = 0.5;

// Note: the anti-kt jet radius parameter R (eta-phi space) is NOT
// redeclared here -- it's read directly from JetParameters.h's `R`
// (included above), the same shared constant JetMatcher.cpp/RecoJets.cpp
// cluster jets with, so this file can never drift out of sync with the
// actual clustering radius used to produce the jetTree it's rendering.

// RGB for jets (truth vs reco) -- distinct from the red/blue hit-energy
// scale and from the cyan EPD-hit color already in use. Alpha (opacity)
// is NOT fixed here -- it's computed per-jet from energy, see
// JET_ALPHA_MIN/MAX and JET_ENERGY_OPACITY_MIN/MAX_GEV below.
static const double RECO_JET_RGB[3] = {1.0, 0.65, 0.0};  // orange
static const double TRUTH_JET_RGB[3] = {0.0, 1.0, 0.0};  // green

// Fallback jet length (mm) for the rare non-forward jet (cos(theta) too
// small/negative to sensibly project onto the front FCS plane).
static const double JET_FALLBACK_LENGTH_MM = 1000.0;

// Minimum jet energy (GeV) to render, matching jet_scale/JetEnergyScaleFineGrid.cpp's
// E0 (its calibration binning starts at 10 GeV; jets below that are never
// included in a bin at all). Not sourced from JetParameters.h since
// JetEnergyScaleFineGrid.cpp doesn't put E0 there either -- it's a
// file-local constant in that script, not a value shared with
// JetMatcher.cpp/RecoJets.cpp the way R is.
static const double MIN_JET_ENERGY_GEV = 10.0;

// Jet opacity scales linearly with energy between these two GeV values
// (clamped outside the range) -- the same [E0, Emax] = [10, 200] GeV
// window jet_scale/JetEnergyScaleFineGrid.cpp bins its calibration map
// over, so the same energy range that already gates which jets render
// at all also sets how opaque they look.
static const double JET_ENERGY_OPACITY_MIN_GEV = MIN_JET_ENERGY_GEV;
static const double JET_ENERGY_OPACITY_MAX_GEV = 200.0;
static const double JET_ALPHA_MIN = 0.3;
static const double JET_ALPHA_MAX = 1.0;

// Track color/opacity: both scale linearly with energy, from
// [0, 0, 1, 0] (0 energy -- fully blue, fully transparent) to
// [1, 0, 0, 1] (this event's max mcparticle energy -- fully red, fully
// opaque). No energy cut and no top-fraction trim on tracks -- every
// mcparticle with a well-defined direction is rendered.

struct Vec3 {
    double x, y, z;
};

// ============================================================
// Energy color scaling (log scale, clamped to [0,1])
// ============================================================

void energy_color(double E, double Emin, double Emax, double out[4]) {
    E = std::max(E, 1e-12);
    Emin = std::max(Emin, 1e-12);
    Emax = std::max(Emax, Emin * 1.0001);

    double fraction = (std::log(E) - std::log(Emin)) / (std::log(Emax) - std::log(Emin));
    fraction = std::min(std::max(fraction, 0.0), 1.0);

    out[0] = fraction;
    out[1] = 0.0;
    out[2] = 1.0 - fraction;
    out[3] = fraction;
}

// ============================================================
// Shift tower from front face to center
// ============================================================

Vec3 shift_block_position(double x, double y, double z, double depth, bool isLeft) {
    double theta = DETECTOR_ANGLE * M_PI / 180.0;
    if (isLeft) theta *= -1;

    double shift = depth / 2.0;
    double dx = std::sin(theta) * shift;
    double dz = std::cos(theta) * shift;

    return { x + dx, y, z + dz };
}

// ============================================================
// Time from origin
// ============================================================

double propagation_time(const Vec3& pos) {
    double distance = std::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
    return distance / C_LIGHT;
}

// ============================================================
// eta -> polar angle theta (radians)
// ============================================================

double eta_to_theta(double eta) {
    return 2.0 * std::atan(std::exp(-eta));
}

// ============================================================
// Static JSON header (identical to FCS2VIRTUE.py's "header" block,
// except tracker_boundary now comes from FCS_Model.json's real
// detector dimensions instead of a placeholder)
// ============================================================

void writeHeader(std::ofstream& f) {
    f << "  \"header\": {\n";
    f << "    \"version\": \"3.1.1\",\n";
    f << "    \"experiment\": \"STAR FCS -- jets: green=truth, orange=reco\",\n";
    f << "    \"energy_unit\": \"GeV\",\n";
    f << "    \"color_bar\": \"Log\",\n";
    f << "    \"scale\": 1.0,\n";
    f << "    \"length_unit\": \"mm\",\n";
    f << "    \"particles\": [\n";
    f << "      {\n";
    f << "        \"size\": 150.0,\n";
    f << "        \"color_rgba\": [1.0, 0.0, 0.0, 1.0],\n";
    f << "        \"ip\": [0.0, 0.0, 0.0],\n";
    f << "        \"angle_rad\": [0.0, 0.0]\n";
    f << "      },\n";
    f << "      {\n";
    f << "        \"size\": 150.0,\n";
    f << "        \"color_rgba\": [1.0, 0.0, 0.0, 1.0],\n";
    f << "        \"ip\": [0.0, 0.0, 0.0],\n";
    f << "        \"angle_rad\": [" << M_PI << ", 0.0]\n";
    f << "      }\n";
    f << "    ],\n";
    f << "    \"tracker_settings\": {\n";
    f << "      \"segment_ns\": " << SEGMENT_NS << ",\n";
    f << "      \"B_field_T\": " << B_FIELD_TESLA << ",\n";
    f << "      \"tracker_boundary\": [" << TRACKER_RADIUS_MM << ", " << TRACKER_Z_MIN_MM << ", " << TRACKER_Z_MAX_MM << "]\n";
    f << "    }\n";
    f << "  },\n";
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_simpletree.root> <n_events> <output.json> [--tracks] [--jets <jettree.root>]" << std::endl;
        return 1;
    }

    const char* input_file = argv[1];
    int n_events_requested = std::atoi(argv[2]);
    const char* output_file = argv[3];

    bool doTracks = false;
    bool doTruthJets = false;
    bool doRecoJets = false;
    std::string jetTreeFile;

    for (int i = 4; i < argc; i++) {
        if (std::strcmp(argv[i], "--tracks") == 0) {
            doTracks = true;
        } else if (std::strcmp(argv[i], "--truth-jets") == 0) {
            doTruthJets = true;
        } else if (std::strcmp(argv[i], "--reco-jets") == 0) {
            doRecoJets = true;
        } else if (std::strcmp(argv[i], "--jets-file") == 0) {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: --jets-file requires a jetTree file path" << std::endl;
                return 1;
            }
            jetTreeFile = argv[++i];
        } else {
            std::cerr << "ERROR: unrecognized option '" << argv[i] << "'" << std::endl;
            return 1;
        }
    }

    if ((doTruthJets || doRecoJets) && jetTreeFile.empty()) {
        std::cerr << "ERROR: --truth-jets/--reco-jets requires --jets-file <jettree.root>" << std::endl;
        return 1;
    }

    // ------------------------------------------------------------
    // SimpleTree (hits, and optionally mcparticle tracks)
    // ------------------------------------------------------------

    TFile f(input_file);
    if (f.IsZombie()) {
        std::cerr << "ERROR: could not open " << input_file << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)f.Get("data");
    if (!tree) {
        std::cerr << "ERROR: no 'data' tree found in " << input_file << std::endl;
        return 1;
    }

    Int_t Cal_nhits;
    Int_t Cal_detid[MAX];
    Float_t Cal_hit_posx[MAX], Cal_hit_posy[MAX], Cal_hit_posz[MAX];
    Float_t Cal_hit_energy[MAX];

    tree->SetBranchAddress("Cal_nhits", &Cal_nhits);
    tree->SetBranchAddress("Cal_detid", Cal_detid);
    tree->SetBranchAddress("Cal_hit_posx", Cal_hit_posx);
    tree->SetBranchAddress("Cal_hit_posy", Cal_hit_posy);
    tree->SetBranchAddress("Cal_hit_posz", Cal_hit_posz);

    bool haveEnergy = (tree->GetBranch("Cal_hit_energy") != nullptr);
    if (haveEnergy) {
        tree->SetBranchAddress("Cal_hit_energy", Cal_hit_energy);
    }

    Int_t mcpart_num;
    Float_t mcpart_px[MAX], mcpart_py[MAX], mcpart_pz[MAX], mcpart_E[MAX];
    Int_t mcpart_charge[MAX];
    Float_t mcpart_Vtx_x[MAX], mcpart_Vtx_y[MAX], mcpart_Vtx_z[MAX];

    if (doTracks) {
        tree->SetBranchAddress("mcpart_num", &mcpart_num);
        tree->SetBranchAddress("mcpart_px", mcpart_px);
        tree->SetBranchAddress("mcpart_py", mcpart_py);
        tree->SetBranchAddress("mcpart_pz", mcpart_pz);
        tree->SetBranchAddress("mcpart_E", mcpart_E);
        tree->SetBranchAddress("mcpart_charge", mcpart_charge);
        tree->SetBranchAddress("mcpart_Vtx_x", mcpart_Vtx_x);
        tree->SetBranchAddress("mcpart_Vtx_y", mcpart_Vtx_y);
        tree->SetBranchAddress("mcpart_Vtx_z", mcpart_Vtx_z);
    }

    // ------------------------------------------------------------
    // jetTree -- truth and reco jets always live in the same file, so
    // this is a single shared TFile/TTree; which branches actually get
    // read depends on doTruthJets/doRecoJets independently.
    // ------------------------------------------------------------

    TFile* jetFile = nullptr;
    TTree* jetTree = nullptr;
    std::vector<float>* truth_pt = nullptr;
    std::vector<float>* truth_eta = nullptr;
    std::vector<float>* truth_phi = nullptr;
    std::vector<float>* truth_E = nullptr;
    std::vector<float>* reco_pt = nullptr;
    std::vector<float>* reco_eta = nullptr;
    std::vector<float>* reco_phi = nullptr;
    std::vector<float>* reco_E = nullptr;

    if (doTruthJets || doRecoJets) {
        jetFile = new TFile(jetTreeFile.c_str());
        if (jetFile->IsZombie()) {
            std::cerr << "ERROR: could not open --jets-file " << jetTreeFile << std::endl;
            return 1;
        }
        jetTree = (TTree*)jetFile->Get("jetTree");
        if (!jetTree) {
            std::cerr << "ERROR: no 'jetTree' found in " << jetTreeFile << std::endl;
            return 1;
        }

        if (doTruthJets) {
            jetTree->SetBranchAddress("truth_pt", &truth_pt);
            jetTree->SetBranchAddress("truth_eta", &truth_eta);
            jetTree->SetBranchAddress("truth_phi", &truth_phi);
            jetTree->SetBranchAddress("truth_E", &truth_E);
        }
        if (doRecoJets) {
            jetTree->SetBranchAddress("reco_pt", &reco_pt);
            jetTree->SetBranchAddress("reco_eta", &reco_eta);
            jetTree->SetBranchAddress("reco_phi", &reco_phi);
            jetTree->SetBranchAddress("reco_E", &reco_E);
        }

        if (jetTree->GetEntries() != tree->GetEntries()) {
            std::cerr << "WARNING: --jets-file jetTree has " << jetTree->GetEntries()
                       << " events but the SimpleTree has " << tree->GetEntries()
                       << " -- are you sure this is the 1:1 corresponding jetTree file?" << std::endl;
        }
    }

    Long64_t total_entries = tree->GetEntries();
    if (jetTree) total_entries = std::min(total_entries, jetTree->GetEntries());
    Long64_t n_events = std::min((Long64_t)n_events_requested, total_entries);
    if (n_events < 0) n_events = 0;

    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "ERROR: could not open output file " << output_file << std::endl;
        return 1;
    }
    out << std::setprecision(6);

    out << "{\n";
    writeHeader(out);
    out << "  \"events\": [\n";

    for (Long64_t event = 0; event < n_events; event++) {
        tree->GetEntry(event);
        if (jetTree) jetTree->GetEntry(event);

        // ----------------------------------------------------
        // First pass: energy range for this event's color scale
        // ----------------------------------------------------
        double Emin = 1.0, Emax = 1.0;
        bool haveRange = false;
        for (int i = 0; i < Cal_nhits; i++) {
            double e = haveEnergy ? Cal_hit_energy[i] : 1.0;
            if (e > 0) {
                if (!haveRange) { Emin = Emax = e; haveRange = true; }
                else { Emin = std::min(Emin, e); Emax = std::max(Emax, e); }
            }
        }

        // ----------------------------------------------------
        // Second pass: build hits (EPD) and blocks (ECAL/HCAL)
        // ----------------------------------------------------
        std::vector<std::string> blockEntries;
        std::vector<std::string> hitEntries;

        for (int i = 0; i < Cal_nhits; i++) {
            int detector = Cal_detid[i];
            double e = haveEnergy ? Cal_hit_energy[i] : 1.0;

            double x = Cal_hit_posx[i] * UNIT_SCALE;
            double y = Cal_hit_posy[i] * UNIT_SCALE;
            double z = Cal_hit_posz[i] * UNIT_SCALE;

            std::ostringstream entry;
            entry << std::setprecision(6);

            if (detector == 0 || detector == 1) {
                // ECAL
                bool isLeft = (detector == 0);
                Vec3 pos = shift_block_position(x, y, z, ECAL_SIZE[2], isLeft);
                double color[4];
                energy_color(e, Emin, Emax, color);
                double eulerY = isLeft ? -DETECTOR_ANGLE : DETECTOR_ANGLE;

                entry << "        {\n";
                entry << "          \"position\": [" << pos.x << ", " << pos.y << ", " << pos.z << "],\n";
                entry << "          \"time_ns\": " << propagation_time(pos) << ",\n";
                entry << "          \"size\": [" << ECAL_SIZE[0] << ", " << ECAL_SIZE[1] << ", " << ECAL_SIZE[2] << "],\n";
                entry << "          \"euler_angles_deg\": [0.0, " << eulerY << ", 0.0],\n";
                entry << "          \"color_rgba\": [" << color[0] << ", " << color[1] << ", " << color[2] << ", " << color[3] << "]\n";
                entry << "        }";
                blockEntries.push_back(entry.str());

            } else if (detector == 2 || detector == 3) {
                // HCAL
                bool isLeft = (detector == 2);
                Vec3 pos = shift_block_position(x, y, z, HCAL_SIZE[2], isLeft);
                double color[4];
                energy_color(e, Emin, Emax, color);
                double eulerY = isLeft ? -DETECTOR_ANGLE : DETECTOR_ANGLE;

                entry << "        {\n";
                entry << "          \"position\": [" << pos.x << ", " << pos.y << ", " << pos.z << "],\n";
                entry << "          \"time_ns\": " << propagation_time(pos) << ",\n";
                entry << "          \"size\": [" << HCAL_SIZE[0] << ", " << HCAL_SIZE[1] << ", " << HCAL_SIZE[2] << "],\n";
                entry << "          \"euler_angles_deg\": [0.0, " << eulerY << ", 0.0],\n";
                entry << "          \"color_rgba\": [" << color[0] << ", " << color[1] << ", " << color[2] << ", " << color[3] << "]\n";
                entry << "        }";
                blockEntries.push_back(entry.str());

            } else if (detector == 4 || detector == 5) {
                // EPD
                Vec3 pos = { x, y, z };

                entry << "        {\n";
                entry << "          \"position\": [" << pos.x << ", " << pos.y << ", " << pos.z << "],\n";
                entry << "          \"time_ns\": " << propagation_time(pos) << ",\n";
                entry << "          \"size\": " << EPD_HIT_SIZE << ",\n";
                entry << "          \"color_rgba\": [0.0, 0.8, 1.0, 0.8]\n";
                entry << "        }";
                hitEntries.push_back(entry.str());
            }
        }

        // ----------------------------------------------------
        // Tracks (mcparticles), optional
        // ----------------------------------------------------
        std::vector<std::string> trackEntries;

        if (doTracks) {
            // Every mcparticle with a well-defined direction is rendered --
            // no energy cut, no top-fraction trim.
            std::vector<int> qualifyingIdx;
            for (int i = 0; i < mcpart_num; i++) {
                double px = mcpart_px[i], py = mcpart_py[i], pz = mcpart_pz[i];
                double p = std::sqrt(px * px + py * py + pz * pz);
                if (p < 1e-9) continue; // at-rest particle, no meaningful direction

                qualifyingIdx.push_back(i);
            }

            // Max energy of ANY mcparticle in this event (not just the
            // rendered/qualifying subset) -- the top of the linear
            // color/opacity scale below.
            double maxE = 0.0;
            for (int i = 0; i < mcpart_num; i++) {
                if (mcpart_E[i] > maxE) maxE = mcpart_E[i];
            }

            for (int i : qualifyingIdx) {
                double px = mcpart_px[i], py = mcpart_py[i], pz = mcpart_pz[i];
                double p = std::sqrt(px * px + py * py + pz * pz);

                double theta = std::acos(pz / p);
                double phi = std::atan2(py, px);
                double qOverP = (mcpart_charge[i] != 0) ? (mcpart_charge[i] / p) : 0.0;

                double vx = mcpart_Vtx_x[i] * UNIT_SCALE;
                double vy = mcpart_Vtx_y[i] * UNIT_SCALE;
                double vz = mcpart_Vtx_z[i] * UNIT_SCALE;

                // Color and opacity both scale linearly with energy, from
                // (blue, transparent) at 0 to (red, opaque) at maxE.
                double trackFrac = (maxE > 0.0) ? (mcpart_E[i] / maxE) : 0.0;
                trackFrac = std::min(std::max(trackFrac, 0.0), 1.0);
                double alpha = trackFrac;

                double rgb[3] = { trackFrac, 0.0, 1.0 - trackFrac };

                // duration_ns: [start, end]. start is the light-travel time
                // from the origin to this track's vertex (same
                // propagation_time() convention used for hits/blocks/jets)
                // -- when the track starts showing up. end is a fixed
                // render-stop time; tracks also stop rendering earlier if
                // they reach the tracker_boundary.
                double startNs = propagation_time({ vx, vy, vz });
                double endNs = 30.0;

                std::ostringstream entry;
                entry << std::setprecision(6);
                entry << "        {\n";
                entry << "          \"qOverP\": " << qOverP << ",\n";
                entry << "          \"angle_rad\": [" << theta << ", " << phi << "],\n";
                entry << "          \"vertex\": [" << vx << ", " << vy << ", " << vz << "],\n";
                entry << "          \"duration_ns\": [" << startNs << ", " << endNs << "],\n";
                entry << "          \"color_rgba\": [" << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << ", " << alpha << "]\n";
                entry << "        }";
                trackEntries.push_back(entry.str());
            }
        }

        // ----------------------------------------------------
        // Jets (truth and/or reco, independently flagged)
        // ----------------------------------------------------
        std::vector<std::string> jetEntries;

        if (doTruthJets || doRecoJets) {
            // vertex = collision point (origin); jets don't carry their
            // own vertex info in jetTree the way tracks do.
            const double jvx = 0.0, jvy = 0.0, jvz = 0.0;

            auto addJets = [&](std::vector<float>* pts, std::vector<float>* etas, std::vector<float>* phis, std::vector<float>* Es, const double rgb[3]) {
                if (!pts) return;
                for (size_t i = 0; i < pts->size(); i++) {
                    double E = (*Es)[i];
                    if (E < MIN_JET_ENERGY_GEV) continue;

                    // Opacity scales linearly with energy across
                    // [JET_ENERGY_OPACITY_MIN_GEV, JET_ENERGY_OPACITY_MAX_GEV],
                    // clamped at both ends.
                    double frac = (E - JET_ENERGY_OPACITY_MIN_GEV) /
                                  (JET_ENERGY_OPACITY_MAX_GEV - JET_ENERGY_OPACITY_MIN_GEV);
                    frac = std::min(std::max(frac, 0.0), 1.0);
                    double alpha = JET_ALPHA_MIN + (JET_ALPHA_MAX - JET_ALPHA_MIN) * frac;

                    double eta = (*etas)[i];
                    double phi = (*phis)[i];
                    double theta = eta_to_theta(eta);

                    // R is defined in eta-phi space, which isn't a fixed
                    // real-space angle: theta = 2*atan(exp(-eta)) gives
                    // dtheta/deta = -sin(theta), and at fixed theta a dphi
                    // displacement sweeps a real angular distance of
                    // sin(theta)*dphi -- so locally, R_etaphi maps to an
                    // actual opening half-angle of R_etaphi * sin(theta)
                    // (smaller for more forward jets, which is the correct
                    // physical behavior: the same eta-phi R is more
                    // collimated in real angle near the beam axis).
                    double R_rad = R * std::sin(theta);

                    double length;
                    double cosTheta = std::cos(theta);
                    if (cosTheta > 0.05) {
                        length = (FCS_FRONT_Z_MM - jvz) / cosTheta;
                    } else {
                        // Not forward-going enough to sensibly reach the
                        // FCS plane (transverse or backward jet) -- draw a
                        // short fixed-length stub instead of a huge/negative
                        // projected length.
                        length = JET_FALLBACK_LENGTH_MM;
                    }

                    // time_ns uses the HALF-length point, not the full
                    // endpoint -- empirically found to look most natural
                    // (the jet "arrives" as light reaches its midpoint,
                    // not its far tip).
                    Vec3 midpoint = {
                        jvx + (length / 2.0) * std::sin(theta) * std::cos(phi),
                        jvy + (length / 2.0) * std::sin(theta) * std::sin(phi),
                        jvz + (length / 2.0) * cosTheta
                    };

                    std::ostringstream entry;
                    entry << std::setprecision(6);
                    entry << "        {\n";
                    entry << "          \"length\": " << length << ",\n";
                    entry << "          \"R_rad\": " << R_rad << ",\n";
                    entry << "          \"angle_rad\": [" << theta << ", " << phi << "],\n";
                    entry << "          \"vertex\": [" << jvx << ", " << jvy << ", " << jvz << "],\n";
                    entry << "          \"time_ns\": " << propagation_time(midpoint) << ",\n";
                    entry << "          \"color_rgba\": [" << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << ", " << alpha << "]\n";
                    entry << "        }";
                    jetEntries.push_back(entry.str());
                }
            };

            if (doTruthJets) addJets(truth_pt, truth_eta, truth_phi, truth_E, TRUTH_JET_RGB);
            if (doRecoJets) addJets(reco_pt, reco_eta, reco_phi, reco_E, RECO_JET_RGB);
        }

        // ----------------------------------------------------
        // Write this event
        // ----------------------------------------------------

        out << "    {\n";
        out << "      \"event_data\": {\n";
        out << "        \"info_text\": \"STAR FCS Event " << event << "\",\n";
        out << "        \"energy_scale\": [" << Emin << ", " << Emax << "]\n";
        out << "      },\n";

        out << "      \"hits\": [\n";
        for (size_t i = 0; i < hitEntries.size(); i++) {
            out << hitEntries[i];
            if (i + 1 < hitEntries.size()) out << ",";
            out << "\n";
        }
        out << "      ],\n";

        out << "      \"blocks\": [\n";
        for (size_t i = 0; i < blockEntries.size(); i++) {
            out << blockEntries[i];
            if (i + 1 < blockEntries.size()) out << ",";
            out << "\n";
        }
        out << "      ]";

        if (doTracks) {
            out << ",\n      \"tracks\": [\n";
            for (size_t i = 0; i < trackEntries.size(); i++) {
                out << trackEntries[i];
                if (i + 1 < trackEntries.size()) out << ",";
                out << "\n";
            }
            out << "      ]";
        }

        if (doTruthJets || doRecoJets) {
            out << ",\n      \"jets\": [\n";
            for (size_t i = 0; i < jetEntries.size(); i++) {
                out << jetEntries[i];
                if (i + 1 < jetEntries.size()) out << ",";
                out << "\n";
            }
            out << "      ]";
        }

        out << "\n    }";
        if (event + 1 < n_events) out << ",";
        out << "\n";
    }

    out << "  ]\n";
    out << "}\n";
    out.close();

    if (jetFile) { jetFile->Close(); delete jetFile; }

    std::cout << "Wrote " << n_events << " events to " << output_file
               << " (tracks=" << (doTracks ? "yes" : "no")
               << ", truth_jets=" << (doTruthJets ? "yes" : "no")
               << ", reco_jets=" << (doRecoJets ? "yes" : "no") << ")" << std::endl;

    return 0;
}
