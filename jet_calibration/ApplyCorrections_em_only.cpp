// Applies an EM-only JES lookup JSON (CalibrationMap_em_only) to JetTrees in
// place (TFile "UPDATE", no undo), adding reco_E_corr and reco_x_F_corr.
//
// scale(x,y,E) = A + B/E + C/E^2 + D exp(-(E-E0)^2 / 2 sigma^2), sigma = sigmaL
// below E0, sigmaR above, evaluated at the jet's own UE-subtracted energy
// E = reco_E - reco_ue_E (numerical inversion): reco_E_corr = E / scale.
// Jets that cannot be put on the calibrated scale get reco_E_corr =
// reco_x_F_corr = -1 and must be dropped by the analysis: no UE estimate
// ("No UE estimate") or outside the map ("Invalid entries": NaN beyond the
// buffered fiducial region). A file without reco_ue_E is an error unless
// "noue" is given (map built without subtraction).

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

#include "JetParameters.h"

struct LookupGrid
{
    int nx, ny;
    std::vector<double> x;
    std::vector<double> y;

    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> C;
    std::vector<std::vector<double>> D;
    std::vector<std::vector<double>> E0;
    std::vector<std::vector<double>> sigmaL;
    std::vector<std::vector<double>> sigmaR;
};

double clamp(double v,double lo,double hi)
{
    if(v<lo) return lo;
    if(v>hi) return hi;
    return v;
}

int find_index(const std::vector<double>& grid,double v)
{
    int N = grid.size();

    if(v <= grid.front()) return 0;
    if(v >= grid.back()) return N-2;

    int lo=0;
    int hi=N-1;

    while(hi-lo>1)
    {
        int mid=(lo+hi)/2;

        if(grid[mid] > v)
            hi = mid;
        else
            lo = mid;
    }

    return lo;
}

double bilinear(
    const std::vector<double>& xg,
    const std::vector<double>& yg,
    const std::vector<std::vector<double>>& Z,
    double x,
    double y
)
{
    int i = find_index(xg,x);
    int j = find_index(yg,y);

    double x1 = xg[i];
    double x2 = xg[i+1];
    double y1 = yg[j];
    double y2 = yg[j+1];

    double q11 = Z[i][j];
    double q12 = Z[i][j+1];
    double q21 = Z[i+1][j];
    double q22 = Z[i+1][j+1];

    double tx = (x - x1)/(x2-x1);
    double ty = (y - y1)/(y2-y1);

    double a = q11*(1-tx) + q21*tx;
    double b = q12*(1-tx) + q22*tx;

    return a*(1-ty) + b*ty;
}

bool invalid(double v)
{
    return !std::isfinite(v);
}

// Minimal JSON parsing (fixed structure)
void read_array_1d(std::ifstream& in, std::vector<double>& out)
{
    out.clear();
    char c;
    double val;

    while(in >> c)
    {
        if(c == '[') break;
    }

    while(in >> val)
    {
        out.push_back(val);
        in >> c;
        if(c == ']') break;
    }
}

void read_array_2d(std::ifstream& in, std::vector<std::vector<double>>& out)
{
    out.clear();
    char c;

    while(in >> c)
    {
        if(c == '[') break;
    }

    while(true)
    {
        in >> c;
        if(c == ']') break;

        if(c == '[')
        {
            std::vector<double> row;
            double val;

            while(in >> val)
            {
                row.push_back(val);
                in >> c;
                if(c == ']') break;
            }

            out.push_back(row);

            in >> c;
            if(c == ']') break;
        }
    }
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " file_list.list calibration.json [noue]" << std::endl;
        return 1;
    }

    std::cout << "Starting ApplyCorrections..." << std::endl;

    gROOT->SetBatch(kTRUE);

    // Floor on the evaluation energy (the model has 1/E terms), well below
    // the 3 GeV bottom of the fit range
    const double E_FLOOR = 1.0;

    std::string list_fname = argv[1];
    std::string calib_json = argv[2];
    const bool subtract_ue = !(argc > 3 && std::string(argv[3]) == "noue");
    long no_ue_counter = 0;
    std::string tree_name = "jetTree";

    std::cout << "Opening JSON: " << calib_json << std::endl;

    std::ifstream in(calib_json.c_str());

    if(!in.is_open())
    {
        std::cerr << "ERROR: Cannot open JSON." << std::endl;
        return 1;
    }

    std::cout << "Reading lookup tables..." << std::endl;

    LookupGrid grid;

    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string json = buffer.str();

    auto extract_int =
    [&](const std::string& key)
    {
        size_t pos = json.find("\"" + key + "\"");
        if(pos == std::string::npos)
        {
            std::cerr << "ERROR: Missing key " << key << std::endl;
            exit(1);
        }

        pos = json.find(':', pos);
        pos++;

        while(pos < json.size() && std::isspace(json[pos]))
            pos++;

        return std::stoi(json.substr(pos));
    };

    // Read nx and ny first
    grid.nx = extract_int("nx");
    grid.ny = extract_int("ny");

    auto extract_array_1d =
    [&](const std::string& key,
        std::vector<double>& out)
    {
        out.clear();

        size_t pos = json.find("\"" + key + "\"");
        if(pos == std::string::npos)
        {
            std::cerr << "ERROR: Missing key " << key << std::endl;
            exit(1);
        }

        pos = json.find('[', pos);
        size_t end = json.find(']', pos);

        std::stringstream ss(
            json.substr(pos + 1, end - pos - 1));

        double value;
        char comma;

        while(ss >> value)
        {
            out.push_back(value);
            ss >> comma;
        }
    };

    extract_array_1d("x_grid", grid.x);
    extract_array_1d("y_grid", grid.y);

    // Sparse per-bin fit params ("bins":[{"i","j","A","B","C","D","E0","sigmaL","sigmaR"},...]);
    // only present over a rectangular (i,j) sub-range, rest left as NaN.
    int imin = grid.nx, imax = -1;
    int jmin = grid.ny, jmax = -1;

    auto init_grid_2d =
    [&](std::vector<std::vector<double>>& v)
    {
        v.assign(grid.nx,
                 std::vector<double>(
                     grid.ny,
                     std::numeric_limits<double>::quiet_NaN()));
    };

    init_grid_2d(grid.A);
    init_grid_2d(grid.B);
    init_grid_2d(grid.C);
    init_grid_2d(grid.D);
    init_grid_2d(grid.E0);
    init_grid_2d(grid.sigmaL);
    init_grid_2d(grid.sigmaR);

    {
        size_t pos = json.find("\"bins\"");
        if(pos == std::string::npos)
        {
            std::cerr << "ERROR: Missing key bins" << std::endl;
            exit(1);
        }

        pos = json.find('[', pos);
        pos++;

        while(true)
        {
            while(pos < json.size() &&
                  (std::isspace(json[pos]) || json[pos] == ','))
                pos++;

            if(pos >= json.size() || json[pos] == ']')
                break;

            if(json[pos] != '{')
            {
                std::cerr << "ERROR: Expected '{' in bins array"
                          << std::endl;
                exit(1);
            }

            size_t obj_start = pos;
            size_t obj_end = json.find('}', obj_start);

            if(obj_end == std::string::npos)
            {
                std::cerr << "ERROR: Unterminated bin object"
                          << std::endl;
                exit(1);
            }

            std::string obj =
                json.substr(obj_start, obj_end - obj_start + 1);

            pos = obj_end + 1;

            auto get_field =
            [&](const std::string& fkey)
            {
                size_t p = obj.find("\"" + fkey + "\"");

                if(p == std::string::npos)
                {
                    std::cerr << "ERROR: bin missing field "
                              << fkey << std::endl;
                    exit(1);
                }

                p = obj.find(':', p);
                p++;

                while(p < obj.size() && std::isspace(obj[p]))
                    p++;

                return std::strtod(obj.c_str() + p, nullptr);
            };

            int bi = (int)get_field("i");
            int bj = (int)get_field("j");

            grid.A[bi][bj]     = get_field("A");
            grid.B[bi][bj]     = get_field("B");
            grid.C[bi][bj]     = get_field("C");
            grid.D[bi][bj]     = get_field("D");
            grid.E0[bi][bj]    = get_field("E0");
            grid.sigmaL[bi][bj] = get_field("sigmaL");
            grid.sigmaR[bi][bj] = get_field("sigmaR");

            imin = std::min(imin, bi);
            imax = std::max(imax, bi);
            jmin = std::min(jmin, bj);
            jmax = std::max(jmax, bj);
        }
    }

    std::cout << "Finished parsing JSON." << std::endl;

    std::cout << "nx = " << grid.nx << std::endl;
    std::cout << "ny = " << grid.ny << std::endl;

    std::cout << "x_grid size = " << grid.x.size() << std::endl;
    std::cout << "y_grid size = " << grid.y.size() << std::endl;

    if(grid.x.empty() || grid.y.empty() || imax < 0 || jmax < 0)
    {
        std::cerr << "ERROR: Lookup tables were not read correctly." << std::endl;
        return 1;
    }

    std::cout << "Valid bin range: i=[" << imin << "," << imax
               << "]  j=[" << jmin << "," << jmax << "]" << std::endl;

    double xmin = grid.x[imin];
    double xmax = grid.x[imax];
    double ymin = grid.y[jmin];
    double ymax = grid.y[jmax];

    std::cout << "Reading file list: " << list_fname << std::endl;

    std::ifstream list_in(list_fname.c_str());
    if(!list_in.is_open())
    {
        std::cerr << "ERROR: Cannot open file list " << list_fname << std::endl;
        return 1;
    }

    std::vector<std::string> root_files;
    std::string fline;
    while(std::getline(list_in, fline))
        if(!fline.empty())
            root_files.push_back(fline);
    list_in.close();

    std::cout << "Found " << root_files.size()
              << " ROOT files." << std::endl;

    long long total_jets = 0;
    long long corrected = 0;
    long long invalid_counter = 0;

    for(size_t fidx = 0; fidx < root_files.size(); fidx++)
    {
        std::cout << "\n==================================================" << std::endl;
        std::cout << "Opening file " << (fidx+1)
                  << "/" << root_files.size() << std::endl;
        std::cout << root_files[fidx] << std::endl;

        TFile f(root_files[fidx].c_str(), "UPDATE");

        if(f.IsZombie())
        {
            std::cout << "Failed to open file." << std::endl;
            continue;
        }

        TTree* tree = (TTree*)f.Get(tree_name.c_str());

        if(!tree)
        {
            std::cout << "Tree not found." << std::endl;
            continue;
        }

        std::cout << "Entries = "
                  << tree->GetEntries() << std::endl;

        std::vector<float>* reco_E = 0;
        std::vector<float>* reco_x = 0;
        std::vector<float>* reco_y = 0;
        std::vector<float>* reco_eta = 0;

        tree->SetBranchAddress("reco_E",&reco_E);
        tree->SetBranchAddress("reco_x",&reco_x);
        tree->SetBranchAddress("reco_y",&reco_y);
        tree->SetBranchAddress("reco_eta",&reco_eta);

        // Per-event sqrt(s) for x_F (fixed simulation value if absent)
        Double_t event_sqrt_s_gev = kSqrtS_MainRun22_GeV;
        if (tree->GetBranch("event_sqrt_s_gev")) tree->SetBranchAddress("event_sqrt_s_gev",&event_sqrt_s_gev);

        std::vector<float>* reco_ue_E = 0;
        std::vector<int>* reco_ue_nproxy = 0;
        if (subtract_ue)
        {
            if (!tree->GetBranch("reco_ue_E") || !tree->GetBranch("reco_ue_nproxy"))
            {
                std::cerr << "ERROR: no reco_ue_E/reco_ue_nproxy branches; run underlying_event/UnderlyingEvent"
                          << " first, or pass noue for a map built without UE subtraction" << std::endl;
                return 1;
            }
            tree->SetBranchAddress("reco_ue_E",&reco_ue_E);
            tree->SetBranchAddress("reco_ue_nproxy",&reco_ue_nproxy);
        }

        std::vector<float> reco_E_corr;
        std::vector<float>* reco_E_corr_ptr = &reco_E_corr;
        TBranch* b_corr =
            tree->Branch("reco_E_corr",&reco_E_corr_ptr);

        // Corrected x_F from the corrected energy and reco_eta
        std::vector<float> reco_x_F_corr;
        std::vector<float>* reco_x_F_corr_ptr = &reco_x_F_corr;
        TBranch* b_xF_corr =
            tree->Branch("reco_x_F_corr",&reco_x_F_corr_ptr);

        Long64_t n_events = tree->GetEntries();

        for(Long64_t ev=0; ev<n_events; ev++)
        {
            if(ev%10000==0)
                std::cout << "Event "
                          << ev << "/"
                          << n_events << std::endl;

            tree->GetEntry(ev);

            if(!reco_E || !reco_x || !reco_y || !reco_eta)
            {
                std::cerr << "ERROR: Null branch pointer at event "
                          << ev << std::endl;
                return 1;
            }

            if(reco_E->size()!=reco_x->size() ||
               reco_E->size()!=reco_y->size() ||
               reco_E->size()!=reco_eta->size())
            {
                std::cerr << "ERROR: Branch sizes differ at event "
                          << ev << std::endl;
                return 1;
            }

            reco_E_corr.clear();
            reco_x_F_corr.clear();

            for(size_t j=0;j<reco_E->size();j++)
            {
                total_jets++;

                double E = reco_E->at(j);
                if (subtract_ue)
                {
                    if (reco_ue_nproxy->at(j) <= 0)
                    {
                        no_ue_counter++;
                        reco_E_corr.push_back(-1.0f);
                        reco_x_F_corr.push_back(-1.0f);
                        continue;
                    }
                    E -= reco_ue_E->at(j);
                }
                // x is reflected positive to match the grid (JetEnergyScaleFineGrid.cpp
                // folds x, not y -- folding y here would look up the wrong bin).
                double x = std::fabs(reco_x->at(j));
                double y = reco_y->at(j);

                x = clamp(x,xmin,xmax);
                y = clamp(y,ymin,ymax);

                double A_grid = bilinear(grid.x,grid.y,grid.A,x,y);
                double B_grid = bilinear(grid.x,grid.y,grid.B,x,y);
                double C_grid = bilinear(grid.x,grid.y,grid.C,x,y);
                double D_grid = bilinear(grid.x,grid.y,grid.D,x,y);
                double E0_grid = bilinear(grid.x,grid.y,grid.E0,x,y);
                double sigmaL_grid = bilinear(grid.x,grid.y,grid.sigmaL,x,y);
                double sigmaR_grid = bilinear(grid.x,grid.y,grid.sigmaR,x,y);

                double E_eval = std::max(E,E_FLOOR);

                double gauss_term = 0.0;

                if(D_grid != 0.0 && std::isfinite(D_grid))
                {
                    double sigma_grid = (E_eval < E0_grid) ? sigmaL_grid : sigmaR_grid;
                    if(sigma_grid != 0.0 && std::isfinite(sigma_grid))
                    {
                        double z = (E_eval - E0_grid)/sigma_grid;
                        gauss_term = D_grid*std::exp(-0.5*z*z);
                    }
                }

                double F_grid =
                    A_grid + B_grid/E_eval + C_grid/(E_eval*E_eval) + gauss_term;

                if(invalid(F_grid) || std::fabs(F_grid)<=1e-8)
                {
                    invalid_counter++;
                    reco_E_corr.push_back(-1.0f);
                    reco_x_F_corr.push_back(-1.0f);
                    continue;
                }

                double E_corr_val = E/F_grid;
                corrected++;
                reco_E_corr.push_back(E_corr_val);
                reco_x_F_corr.push_back(
                    (float)computeFeynmanX(E_corr_val, reco_eta->at(j), event_sqrt_s_gev));
            }

            b_corr->Fill();
            b_xF_corr->Fill();
        }

        std::cout << "Writing tree..." << std::endl;
        tree->Write("", TObject::kOverwrite);

        std::cout << "Closing file..." << std::endl;
        f.Close();
    }

    std::cout << "\nFinished." << std::endl;
    std::cout << "Total jets      : " << total_jets << std::endl;
    std::cout << "Corrected       : " << corrected << std::endl;
    std::cout << "Invalid entries : " << invalid_counter << std::endl;
    if (subtract_ue)
        std::cout << "No UE estimate  : " << no_ue_counter << std::endl;

    return 0;
}