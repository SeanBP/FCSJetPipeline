// Hot / dead FCS channel mask per run period (Sean 2026-11-24), shared by
// StFcsMuToEventMaker (masked hits get zero energy before reclustering, so
// they enter no cluster or point) and StSimpleReaderMaker (masked hits are left
// out of the hit list and the offline trigger 2x2 sums).
//
// The file is asymmetry/channel_qa/channel_mask.cpp's channel_mask.txt (or any
// file with the same leading columns): one row per masked channel and period,
//   period first_run last_run det id ...
// (# comments; columns after id are ignored). A run is masked for a channel if
// first_run <= run <= last_run on one of its rows. No file = nothing masked.
#ifndef FcsChannelMask_h
#define FcsChannelMask_h

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class FcsChannelMask
{
 public:
  FcsChannelMask() : mCachedRun(-1) {}
  // Reads the mask file; returns the number of rows read (-1 if it cannot be opened)
  int Load(const char* file)
  {
    mRanges.clear();
    std::ifstream in(file);
    if (!in) return -1;
    std::string line;
    int n = 0;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream s(line);
      int period, first, last, det, id;
      if (!(s >> period >> first >> last >> det >> id)) continue;
      mRanges[std::make_pair(det, id)].push_back(std::make_pair(first, last));
      n++;
    }
    mCachedRun = -1;
    return n;
  }
  bool Empty() const { return mRanges.empty(); }
  // Is (det, id) masked in run?
  bool Masked(int run, int det, int id)
  {
    if (mRanges.empty()) return false;
    if (run != mCachedRun)
    {
      mCurrent.clear();
      for (std::map<std::pair<int, int>, std::vector<std::pair<int, int> > >::const_iterator it = mRanges.begin(); it != mRanges.end(); ++it)
        for (size_t k = 0; k < it->second.size(); k++)
          if (it->second[k].first <= run && run <= it->second[k].second) { mCurrent[it->first] = true; break; }
      mCachedRun = run;
    }
    return mCurrent.count(std::make_pair(det, id)) > 0;
  }
  // Masked channels in run (after a Masked() call for that run)
  size_t NumMaskedInRun() const { return mCurrent.size(); }

 private:
  std::map<std::pair<int, int>, std::vector<std::pair<int, int> > > mRanges;
  std::map<std::pair<int, int>, bool> mCurrent;
  int mCachedRun;
};

#endif
