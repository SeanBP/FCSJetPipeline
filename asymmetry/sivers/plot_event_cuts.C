// Event-level cut variables with their cuts (Sean 2026-11-17): for each event
// cut of AnHistograms, the distribution of the variable cut on (N-1: every
// other event cut applied), a solid vertical line at the nominal cut and
// dotted lines at its systematic variations, and the share of events each removes.
//   qa_cut_vz            measured vertex z: |vz| < 100 cm (nominal), < 60 cm (variation)
//   qa_cut_vz_source     vertex source: the "measured vertex only" variation
//                        drops the events at the nominal z (no measurement)
//   qa_cut_bbc_e_ntiles  east BBC tiles fired (FMS definition): >= 1
//   qa_cut_tof_mult      TOF trigger multiplicity: >= 3 (> 2)
// The cut positions come from the TNamed "event_cut_lines" AnHistograms writes
// (AnConfig.h / JetParameters.h values), so the lines follow the code.
// Input: one yield file, an AnMerge output, or a .list of yield files (summed).
// Writes <prefix>_<variable>.png and <prefix>_all.png (2x2), and a summary
// <prefix>.txt (events and the share removed by each cut and variation).
// Usage: root -b -q 'plot_event_cuts.C+("merged.root","event_cuts")'
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

struct CutLines { std::vector<double> nom, var; TString varlabel; };

static std::vector<double> ParseEdges(const TString &s)
{
    std::vector<double> v;
    TObjArray *a = s.Tokenize(",");
    for (int i = 0; i < a->GetEntriesFast(); i++) v.push_back(((TObjString *)a->At(i))->GetString().Atof());
    delete a;
    return v;
}

void plot_event_cuts(const char *input, const char *prefix = "event_cuts")
{
    gStyle->SetOptStat(0);
    const char *names[4] = {"qa_cut_vz", "qa_cut_vz_source", "qa_cut_bbc_e_ntiles", "qa_cut_tof_mult"};
    TH1D *h[4] = {0, 0, 0, 0};
    TString lines_str;
    // one file or a list of files
    std::vector<std::string> files;
    if (TString(input).EndsWith(".list"))
    {
        std::ifstream in(input);
        std::string l;
        while (in >> l) if (!l.empty() && l[0] != '#') files.push_back(l);
    }
    else files.push_back(input);
    int nread = 0;
    for (size_t i = 0; i < files.size(); i++)
    {
        TFile *f = TFile::Open(files[i].c_str());
        if (!f || f->IsZombie()) { delete f; continue; }
        bool any = false;
        for (int k = 0; k < 4; k++)
        {
            TH1D *x = (TH1D *)f->Get(names[k]);
            if (!x) continue;
            any = true;
            if (!h[k]) { h[k] = (TH1D *)x->Clone(); h[k]->SetDirectory(0); }
            else h[k]->Add(x);
        }
        if (lines_str == "") if (TNamed *n = (TNamed *)f->Get("event_cut_lines")) lines_str = n->GetTitle();
        if (any) nread++;
        f->Close();
        delete f;
    }
    if (!nread) { printf("plot_event_cuts: no qa_cut_* histograms in %s (yield files from 2026-11-17 on)\n", input); return; }
    // cut positions: "<hist>:<nominal>:<variations>[:<label>];..."
    std::map<std::string, CutLines> cuts;
    if (lines_str == "") printf("plot_event_cuts: WARNING no event_cut_lines in the input; no cut lines drawn\n");
    TObjArray *ent = lines_str.Tokenize(";");
    for (int i = 0; i < ent->GetEntriesFast(); i++)
    {
        TString e = ((TObjString *)ent->At(i))->GetString();
        TObjArray *p = e.Tokenize(":");
        // Tokenize drops empty fields: read them by position from the string instead
        delete p;
        std::vector<TString> fld;
        Ssiz_t st = 0;
        for (;;)
        {
            Ssiz_t c = e.Index(":", st);
            if (c == kNPOS) { fld.push_back(e(st, e.Length() - st)); break; }
            fld.push_back(e(st, c - st));
            st = c + 1;
        }
        CutLines cl;
        if (fld.size() > 1) cl.nom = ParseEdges(fld[1]);
        if (fld.size() > 2) cl.var = ParseEdges(fld[2]);
        if (fld.size() > 3) cl.varlabel = fld[3];
        cuts[fld[0].Data()] = cl;
    }
    delete ent;

    FILE *txt = fopen(Form("%s.txt", prefix), "w");
    fprintf(txt, "# event-level cut variables (plot_event_cuts.C), input %s (%d files); N-1 distributions\n", input, nread);
    fprintf(txt, "# variable events_in_plot cut removed_share [variation removed_share]\n");
    // share of the plotted events outside [lo, hi) for a set of edges (2 edges = window, 1 = lower edge)
    auto removed = [](TH1D *x, const std::vector<double> &e)
    {
        const double tot = x->Integral(0, x->GetNbinsX() + 1);
        if (!(tot > 0) || e.empty()) return 0.0;
        double kept = 0;
        for (int b = 0; b <= x->GetNbinsX() + 1; b++)
        {
            const double c = x->GetBinCenter(b);
            const bool in = e.size() >= 2 ? (c > e[0] && c < e[1]) : c > e[0];
            if (in) kept += x->GetBinContent(b);
        }
        return 1 - kept / tot;
    };
    const char *label[4] = {"vertex z (measured)", "vertex source", "east BBC tiles", "TOF multiplicity"};
    bool write_txt = true;  // the summary once (the 2x2 page redraws)
    auto draw_one = [&](int k, TVirtualPad *pad)
    {
        pad->cd();
        pad->SetLeftMargin(0.13); pad->SetBottomMargin(0.12); pad->SetTopMargin(0.12);
        TH1D *x = h[k];
        if (!x) { TLatex t; t.SetNDC(); t.DrawLatex(0.2, 0.5, Form("%s: not in the input", names[k])); return; }
        x->SetLineColor(kBlack); x->SetLineWidth(2); x->SetFillColor(kGray); x->SetFillStyle(1001);
        x->SetTitle("");
        x->GetXaxis()->SetTitleSize(0.045); x->GetYaxis()->SetTitleSize(0.045); x->GetYaxis()->SetTitleOffset(1.3);
        const bool cat = k == 1;
        if (!cat) pad->SetLogy();
        x->SetMinimum(cat ? 0 : 0.5);
        const double hmax = x->GetBinContent(x->GetMaximumBin());  // not GetMaximum(): set on the first draw
        x->SetMaximum(cat ? 1.8 * hmax : 40 * hmax);
        x->Draw("hist");
        // lines stop below the legend
        const double ylo = cat ? 0 : 0.5, yhi = cat ? 1.1 * hmax : 2.5 * hmax;
        TLegend *lg = new TLegend(0.42, 0.64, 0.89, 0.80);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
        lg->AddEntry(x, Form("events (%.0f)", x->Integral(0, x->GetNbinsX() + 1)), "f");
        const double tot = x->Integral(0, x->GetNbinsX() + 1);
        if (cat)
        {
            // the measured-only variation keeps EPD and BBC TAC (bins 2, 3)
            TLine *l = new TLine(3, ylo, 3, yhi); l->SetLineColor(kBlue + 1); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
            const double rm = tot > 0 ? (x->GetBinContent(1) + x->GetBinContent(4)) / tot : 0;
            lg->AddEntry(l, Form("measured vertex only: removes %.1f%%", 100 * rm), "l");
            lg->AddEntry((TObject *)0, "nominal: all sources kept", "");
            if (write_txt) fprintf(txt, "vz_source %.0f none 0 measured_only %.4f\n", tot, rm);
        }
        else
        {
            const CutLines &cl = cuts[names[k]];
            const double rn = removed(x, cl.nom);
            for (size_t i = 0; i < cl.nom.size(); i++)
            {
                TLine *l = new TLine(cl.nom[i], ylo, cl.nom[i], yhi); l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->Draw();
                if (i == 0) lg->AddEntry(l, Form("cut: removes %.2f%%", 100 * rn), "l");
            }
            const double rv = removed(x, cl.var);
            for (size_t i = 0; i < cl.var.size(); i++)
            {
                TLine *l = new TLine(cl.var[i], ylo, cl.var[i], yhi); l->SetLineColor(kBlue + 1); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
                if (i == 0) lg->AddEntry(l, Form("variation%s%s: removes %.2f%%", cl.varlabel != "" ? " " : "", cl.varlabel.Data(), 100 * rv), "l");
            }
            if (write_txt)
            {
                fprintf(txt, "%s %.0f %s %.4f", names[k] + 7, tot, cl.nom.empty() ? "none" : "nominal", rn);
                if (!cl.var.empty()) fprintf(txt, " variation %.4f", rv);
                fprintf(txt, "\n");
            }
        }
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.04);
        t.DrawLatex(0.13, 0.92, Form("Run 22 FCS: %s", label[k]));
    };
    // keep the N-1 condition text (the histogram title) before blanking it
    TString cond[4];
    for (int k = 0; k < 4; k++) if (h[k]) cond[k] = h[k]->GetTitle();
    auto note = [&](int k)
    {
        TLatex t; t.SetNDC(); t.SetTextSize(0.03); t.SetTextColor(kGray + 2);
        t.DrawLatex(0.15, 0.835, "N-1: " + cond[k]);
    };
    const char *suf[4] = {"vz", "vz_source", "bbc_e_ntiles", "tof_mult"};
    for (int k = 0; k < 4; k++)
    {
        TCanvas c(Form("c%d", k), "", 900, 700);
        draw_one(k, &c);
        note(k);
        c.SaveAs(Form("%s_%s.png", prefix, suf[k]));
    }
    write_txt = false;
    TCanvas ca("call", "", 1600, 1200);
    ca.Divide(2, 2);
    for (int k = 0; k < 4; k++) { draw_one(k, ca.cd(k + 1)); note(k); }
    ca.SaveAs(Form("%s_all.png", prefix));
    fclose(txt);
    printf("plot_event_cuts: %d files -> %s_{vz,vz_source,bbc_e_ntiles,tof_mult,all}.png, %s.txt\n", nread, prefix, prefix);
}
