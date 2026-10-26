#ifndef PROXY_SEARCH_H
#define PROXY_SEARCH_H

// The proxy-event search: for a jet, average the ECal point energy inside a
// cone at the jet's position over other events that pass ProxySelection.h.
// Shared by UnderlyingEvent.cpp and UnderlyingEventPlayground/.
//
// Add every event to a ProxyPool in time order with its group (run in data,
// ptHat bin in simulation), highest EM trigger level, reco jets and ECal
// points, then call Estimate() per jet.
//
// Jets are indexed by group, EM level and 5%-wide pT bin. Candidates from the
// bins overlapping the pT window are visited nearest in event number first. A
// candidate event is accepted if it is a different, not yet used event; its
// jet passes the pT window, fiducial region, azimuthal band and separation;
// and none of its reco jets lies within the separation of the UE cone (unless
// built with PROXY_JET_VETO=0, ProxySelection.h).

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "JetParameters.h"
#include "ProxySelection.h"

struct ProxyJet { float pt, eta, phi, x, y; };

// Band half-width against UE position, tabulated once on a grid over one half
// of the reco fiducial region (the detector is left-right symmetric)
struct ProxyWidthTable
{
    double x_lo, y_lo, step; int nx, ny;
    std::vector<float> width;
    void Build(double grid_cm)
    {
        const FiducialRect &r = kFiducialRectEcal;
        const double b = reco_fiducial_buffer_ecal;
        step = grid_cm; x_lo = r.x_inner + b; y_lo = r.y_min + b;
        nx = (int)std::ceil((r.x_outer - b - x_lo) / step); ny = (int)std::ceil((r.y_max - b - y_lo) / step);
        width.assign(nx * ny, 0.f);
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            {
                double eta, phi;
                computeEtaPhi(x_lo + (i + 0.5) * step, y_lo + (j + 0.5) * step, eta, phi);
                width[i * ny + j] = ProxyHalfWidthDeg(eta, phi);
            }
    }
    double Get(double x, double y) const
    {
        const int i = std::min(nx - 1, std::max(0, (int)((std::fabs(x) - x_lo) / step)));
        const int j = std::min(ny - 1, std::max(0, (int)((y - y_lo) / step)));
        return width[i * ny + j];
    }
};

struct ProxyResult
{
    int n;                          // proxy events used; -1 if the jet cannot be estimated at all
    int n_nominal;                  // how many of them came from the nominal pT window
    double window;                  // pT window half-width finally used, as a fraction of the jet pT
    double ue_E, ue_pt;             // mean cone energy and mean cone sum of point pT; -1 if n <= 0
    long n_scanned;                 // candidates examined
    int reach;                      // distance in pool events to the farthest proxy event
    std::vector<double> E_list;     // cone energy of each proxy event, in the order found
    double jet_eta, jet_dr, jet_dphi; // proxy jets' mean eta, distance to the cone and |dphi| to it (deg); -1 if n <= 0
    std::vector<int> proxy_event, proxy_jet; // pool indices of each proxy event and its proxy jet (GetEvent, GetJet), as E_list
    double half_width;              // band half-width used (deg), from the jet's position; -1 if not estimated
    ProxyResult() : n(-1), n_nominal(0), window(proxy_pt_window_fraction), ue_E(-1.0), ue_pt(-1.0), n_scanned(0), reach(0),
                    jet_eta(-1.0), jet_dr(-1.0), jet_dphi(-1.0), half_width(-1.0) {}
};

class ProxyPool
{
public:
    struct Event
    {
        int group, em_level;
        int jet_first, jet_n;       // range in jets_
        int pt_first, pt_n;         // range in the point arrays
    };

    explicit ProxyPool(double table_grid_cm = 2.0) : require_em_trigger_(true), jet_veto_(proxy_jet_veto), pt_min_(reco_jet_pt_min_em) { table_.Build(table_grid_cm); }

    // Minimum-bias mode (samples without trigger or jet pT cuts, e.g. the JES
    // calibration sample): every fiducial jet is estimated and can be a proxy
    // jet, and untriggered events match untriggered ones. Call before AddEvent.
    void SetMinimumBias() { require_em_trigger_ = false; pt_min_ = 0.0; }

    // The jet veto for this pool; proxy_jet_veto (PROXY_JET_VETO) unless set
    // (FMSJetPipeline's FMS test runs pools with and without it)
    void SetJetVeto(bool on) { jet_veto_ = on; }

    // Adds one event and returns its index. Events must be added in time
    // order. Points are ECal points that already passed the energy cut.
    int AddEvent(int group, int em_level, const std::vector<ProxyJet> &jets,
                 const std::vector<float> &point_eta, const std::vector<float> &point_phi, const std::vector<float> &point_E)
    {
        Event e;
        e.group = group; e.em_level = em_level;
        e.jet_first = (int)jets_.size(); e.jet_n = (int)jets.size();
        e.pt_first = (int)p_E_.size(); e.pt_n = (int)point_E.size();
        const int index = (int)events_.size();
        for (size_t j = 0; j < jets.size(); j++)
        {
            if ((em_level >= 0 || !require_em_trigger_) && jets[j].pt >= pt_min_ && jets[j].pt > 0)
            {
                Ref r = {index, (int)jets_.size()};
                index_[Key(group, em_level, PtBin(jets[j].pt))].push_back(r);
            }
            jets_.push_back(jets[j]);
        }
        p_eta_.insert(p_eta_.end(), point_eta.begin(), point_eta.end());
        p_phi_.insert(p_phi_.end(), point_phi.begin(), point_phi.end());
        p_E_.insert(p_E_.end(), point_E.begin(), point_E.end());
        events_.push_back(e);
        return index;
    }

    int NEvents() const { return (int)events_.size(); }
    size_t NJets() const { return jets_.size(); }
    size_t NPoints() const { return p_E_.size(); }
    size_t NIndexBins() const { return index_.size(); }
    const Event &GetEvent(int i) const { return events_[i]; }
    const ProxyJet &GetJet(int i) const { return jets_[i]; }
    double MemoryMB() const
    {
        return (events_.size() * sizeof(Event) + jets_.size() * (sizeof(ProxyJet) + sizeof(Ref)) + 3.0 * p_E_.size() * sizeof(float)) / 1048576.0;
    }

    // UE estimate for a jet. n = -1 if the event fired no EM trigger, the jet is
    // below the reco pT threshold, or it is outside the fiducial region (in
    // minimum-bias mode only the last applies). n_max <= 0: no limit.
    ProxyResult Estimate(int event, int jet_in_event, int n_max = proxy_n_max) const
    {
        const int n_target = n_max > 0 ? n_max : 2000000000;
        ProxyResult res;
        const Event &ev = events_[event];
        const ProxyJet &jet = jets_[ev.jet_first + jet_in_event];
        if ((ev.em_level < 0 && require_em_trigger_) || jet.pt < pt_min_ || !(jet.pt > 0) || !ProxyInFiducial(jet.x, jet.y))
            return res;
        res.n = 0;

        const double half_width = table_.Get(jet.x, jet.y);
        res.half_width = half_width;
        std::vector<int> used_events;
        std::vector<double> pt_list;
        double sum_eta = 0.0, sum_dr = 0.0, sum_dphi = 0.0;

        // The nominal window, up to the limit if any; then, only for a jet left
        // below the minimum, wider windows until the minimum is reached
        const int n_min = std::min(proxy_n_min, n_target);
        Search(event, jet, half_width, res.window, n_target, res, used_events, pt_list, sum_eta, sum_dr, sum_dphi);
        res.n_nominal = (int)res.E_list.size();
        while ((int)res.E_list.size() < n_min && res.window + proxy_pt_window_step <= proxy_pt_window_max + 1e-9)
        {
            res.window += proxy_pt_window_step;
            Search(event, jet, half_width, res.window, n_min, res, used_events, pt_list, sum_eta, sum_dr, sum_dphi);
        }

        res.n = (int)res.E_list.size();
        if (res.n > 0)
        {
            res.ue_E = 0.0; res.ue_pt = 0.0;
            for (int i = 0; i < res.n; i++) { res.ue_E += res.E_list[i]; res.ue_pt += pt_list[i]; }
            res.ue_E /= res.n; res.ue_pt /= res.n;
            res.jet_eta = sum_eta / res.n; res.jet_dr = sum_dr / res.n; res.jet_dphi = sum_dphi / res.n;
        }
        return res;
    }

private:
    struct Ref { int event, jet; };              // one indexed jet
    static bool ByEvent(const Ref &a, const Ref &b) { return a.event < b.event; }
    typedef std::pair<long long, int> KeyType;
    static KeyType Key(int group, int em_level, int pt_bin) { return KeyType((long long)group * 4 + em_level, pt_bin); }
    // 5% wide; bin 0 starts at reco_jet_pt_min_em, lower pT gives negative bins
    static int PtBin(double pt) { return (int)std::floor(std::log(pt / reco_jet_pt_min_em) / std::log(1.05)); }
    bool require_em_trigger_, jet_veto_;
    double pt_min_;

    // Adds proxy events whose proxy jet is inside the given pT window, nearest
    // in event number first, until there are `want` in all or the group is exhausted
    void Search(int event, const ProxyJet &jet, double half_width, double window, int want,
                ProxyResult &res, std::vector<int> &used_events, std::vector<double> &pt_list,
                double &sum_eta, double &sum_dr, double &sum_dphi) const
    {
        const Event &ev = events_[event];

        // The index bins overlapping the pT window, each with a cursor either side of this event
        std::vector<const std::vector<Ref> *> bins;
        std::vector<int> left, right;
        for (int b = PtBin(std::max(jet.pt * (1.0 - window), std::max(pt_min_, 1e-3))); b <= PtBin(jet.pt * (1.0 + window)); b++)
        {
            std::map<KeyType, std::vector<Ref> >::const_iterator it = index_.find(Key(ev.group, ev.em_level, b));
            if (it == index_.end()) continue;
            const Ref key = {event, 0};
            const int pos = (int)(std::lower_bound(it->second.begin(), it->second.end(), key, ByEvent) - it->second.begin());
            bins.push_back(&it->second); left.push_back(pos - 1); right.push_back(pos);
        }

        while ((int)res.E_list.size() < want)
        {
            // Nearest unvisited candidate in event number, over all bins and both directions
            int best_bin = -1, best_dist = 0; bool best_right = false;
            for (size_t b = 0; b < bins.size(); b++)
            {
                if (left[b] >= 0)
                {
                    const int d = event - (*bins[b])[left[b]].event;
                    if (best_bin < 0 || d < best_dist) { best_bin = (int)b; best_dist = d; best_right = false; }
                }
                if (right[b] < (int)bins[b]->size())
                {
                    const int d = (*bins[b])[right[b]].event - event;
                    if (best_bin < 0 || d < best_dist) { best_bin = (int)b; best_dist = d; best_right = true; }
                }
            }
            if (best_bin < 0) break; // the group is exhausted
            const Ref cand = best_right ? (*bins[best_bin])[right[best_bin]++] : (*bins[best_bin])[left[best_bin]--];
            res.n_scanned++;

            if (cand.event == event) continue;
            const ProxyJet &pj = jets_[cand.jet];
            if (!ProxyPtAllowed(jet.pt, pj.pt, window) || !ProxyInFiducial(pj.x, pj.y)) continue;
            if (!ProxyJetAllowed(jet.eta, jet.phi, half_width, pj.eta, pj.phi)) continue;
            if (std::find(used_events.begin(), used_events.end(), cand.event) != used_events.end()) continue;

            // Only jets above the analysis threshold veto (as in data, where
            // softer jets do not exist); in minimum-bias mode they are UE.
            const Event &pe = events_[cand.event];
            bool clear = true;
            for (int o = pe.jet_first; o < pe.jet_first + pe.jet_n && clear && jet_veto_; o++)
                if (jets_[o].pt >= reco_jet_pt_min_em &&
                    ProxyDeltaR(jets_[o].eta, jets_[o].phi, jet.eta, jet.phi) < proxy_min_dr) clear = false;
            if (!clear) continue;

            double e_ue = 0.0, pt_ue = 0.0;
            for (int i = pe.pt_first; i < pe.pt_first + pe.pt_n; i++)
                if (ProxyDeltaR(p_eta_[i], p_phi_[i], jet.eta, jet.phi) < R)
                {
                    e_ue += p_E_[i];
                    pt_ue += p_E_[i] / std::cosh(p_eta_[i]);
                }
            res.E_list.push_back(e_ue);
            res.proxy_event.push_back(cand.event);
            res.proxy_jet.push_back(cand.jet);
            pt_list.push_back(pt_ue);
            used_events.push_back(cand.event);
            sum_eta += pj.eta;
            sum_dr += ProxyDeltaR(pj.eta, pj.phi, jet.eta, jet.phi);
            sum_dphi += std::fabs(std::remainder(pj.phi - jet.phi, 2.0 * M_PI)) * 180.0 / M_PI;
            res.reach = std::max(res.reach, best_dist);
        }
    }

    ProxyWidthTable table_;
    std::vector<Event> events_;
    std::vector<ProxyJet> jets_;
    std::vector<float> p_eta_, p_phi_, p_E_;
    std::map<KeyType, std::vector<Ref> > index_;
};

#endif
