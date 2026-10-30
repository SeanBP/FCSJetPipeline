// Dumps the conditions-DB ECal gain and gain correction of every tower (1496)
// at a run's start time, in both of X. Liang's calibration file formats, to
// compare with his pi0 anchor files:
//   db_fcsgaincorr_<run>.txt         "ehp ns dep ch gaincorr" (as fcsgaincorr_<day>_<N>.txt)
//   db_fcsgaincorrOffline_<run>.txt  "ns id gaincorr"          (as fcsgaincorrOffline_<day>_<N>.txt)
//   db_fcsgain_<run>.txt             "ns id gain gaincorr gainOnline etGain"
// run = 0: no run-based date; date/time (yyyymmdd, hhmmss) pinned if given,
// e.g. 20230410 0, the "y2023" startup date the simulation chain uses
// (sim_to_jet runSimBfc.C with calibRun = 0; output files then named
// <run>_<date>).
// The DB date comes from GainPlayground/get_fcs_gain.C's RunNumberToDateTime().
//
// Usage (STAR environment): root4star -b -q 'dump_ecal_gaincorr.C(23007048)'
void dump_ecal_gaincorr(int run = 23007048, int date = 0, int time = 0)
{
    gROOT->LoadMacro("../../GainPlayground/get_fcs_gain.C"); // RunNumberToDateTime(), and the chain pointer
    gROOT->LoadMacro("bfc.C");
    bfc(-1, "y2023,db,fcsDb", "");
    St_db_Maker *dbMk = (St_db_Maker *)chain->GetMaker("db");
    int idat = 0, itim = 0;
    if (run != 0)
    {
        RunNumberToDateTime(run, idat, itim);
        dbMk->SetDateTime(idat, itim);
    }
    else if (date != 0)
    {
        idat = date; itim = time;
        dbMk->SetDateTime(idat, itim);
    }
    TString tag = Form("%d", run);
    if (run == 0) tag = Form("0_%d", idat);
    chain->Init();
    StFcsDbMaker *fcsdbmkr = (StFcsDbMaker *)chain->GetMaker("fcsDbMkr");
    fcsdbmkr->InitRun(run != 0 ? run : 1);
    StFcsDb *db = (StFcsDb *)chain->GetDataSet("fcsDb");

    ofstream fdc(Form("db_fcsgaincorr_%s.txt", tag.Data())), foff(Form("db_fcsgaincorrOffline_%s.txt", tag.Data())), fg(Form("db_fcsgain_%s.txt", tag.Data()));
    for (int det = 0; det < 2; det++)
        for (int id = 0; id < db->maxId(det); id++)
        {
            int ehp, ns, crt, slt, dep, ch;
            db->getDepfromId(det, id, ehp, ns, crt, slt, dep, ch);
            const float g = db->getGain(det, id), c = db->getGainCorrection(det, id);
            fdc << " " << ehp << "  " << ns << "    " << dep << "    " << ch << "     " << c << endl;
            foff << ns << "    " << id << "     " << c << endl;
            fg << ns << " " << id << " " << g << " " << c << " " << db->getGainOnline(det, id) << " " << db->getEtGain(det, id, 1.0) << endl;
        }
    printf("dump_ecal_gaincorr: run %d, DB date %d %06d, wrote db_fcsgain(corr)(Offline)_%s.txt\n", run, idat, itim, tag.Data());
}
