// Keep only the vertex branches of a SimpleTree (vertex test)
void slim(const char *in, const char *out)
{
    TFile f(in); TTree *t = (TTree *)f.Get("data"); if (!t) { f.ls(); return; }
    t->SetBranchStatus("*", 0);
    const char *b[] = {"event_run","event_id","Spin_config","event_vertex_z","event_vertex_source","bbc_tac_diff","bbc_vz_tac","bbc_vz_det","epd_vz","vpd_vz","epd_avgtac_e","epd_avgtac_w","epd_ntac_e","epd_ntac_w","Trig_flag","fcs_tcu_bits"};
    for (int i = 0; i < 16; i++) if (t->GetBranch(b[i])) t->SetBranchStatus(b[i], 1);
    TFile o(out, "RECREATE"); TTree *s = t->CloneTree(-1, "fast"); s->Write(); o.Close();
}
