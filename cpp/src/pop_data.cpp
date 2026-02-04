#include "pop_data.h"
#include "stats.h"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip>   // add near the top of pop_data.cpp if not already included
#include <ostream>










PopData::PopData(const std::string& path, const std::string& name, float tailSize, int regRange) :
  name_(name), path_(path), tailSize_(tailSize), regRange_(regRange) 
  {
  std::ifstream in(path_);
  if (!in) throw std::runtime_error("cannot open pop file");

  v.clear();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    double prs, pheno;
    if (!(iss >> prs >> pheno)) continue;   // skips header + any junk line

    v.push_back({prs, pheno, -1, -1, 0.0});
  }

  if (v.empty())
    throw std::runtime_error("POP file read zero numeric rows");

  // Global PRS normalization (match Python: np.std ddof=0)
  double m = 0.0;
  for (const auto& p : v) m += p.prs;
  m /= (double)v.size();

  double var = 0.0;
  for (const auto& p : v) {
    double d = p.prs - m;
    var += d * d;
  }
  var /= (double)v.size(); // ddof=0

  double pStd = std::sqrt(var);
  if (pStd == 0.0)
    throw std::runtime_error("PRS standard deviation is zero");

  for (auto& p : v) p.prs /= pStd;
  sort_and_rank();



  double sumx=0, sumy=0, sumxx=0, sumyy=0, sumxy=0;
  const double n = (double)v.size();

  for (const auto& p : v) {
    const double x = p.prs;
    const double y = p.pheno;
    sumx  += x;
    sumy  += y;
    sumxx += x*x;
    sumyy += y*y;
    sumxy += x*y;
	}

  const double r = corr_sums(n, sumx, sumy, sumxx, sumyy, sumxy);
  r2_ = std::isfinite(r) ? r*r : std::numeric_limits<double>::quiet_NaN();

  
  
  
}











void PopData::sort_and_rank() {
  std::sort(v.begin(), v.end(),
    [](const PopPoint& x, const PopPoint& y) {
      return x.pheno < y.pheno;
    });

  int n = v.size();
  for (int i = 0; i < n; i++) {
    v[i].rank = i;
    int p = int(100.0 * i / n);
    if (p > 99) p = 99;
    v[i].pct = p;
  }
}

void PopData::run_reg(int width) {
  int lo = (100 - width) / 2;
  int hi = 99 - lo;

  std::vector<double> x, y;
  for (auto& p : v)
    if (p.pct >= lo && p.pct <= hi) {
      x.push_back(p.pheno);
      y.push_back(p.prs);
    }

  reg_with_intercept(x, y, a, b);

  for (auto& p : v)
    p.expected = a + b * p.pheno;
}

void PopData::calculate_means() {
  bins.assign(100, BinStats());

  for (auto& p : v) {
    auto& b = bins[p.pct];
    b.n++;
    b.mean_prs += p.prs;
    b.mean_exp += p.expected;
  }

  for (auto& b : bins)
    if (b.n > 0) {
      b.mean_prs /= b.n;
      b.mean_exp /= b.n;
    }

  for (int pct = 0; pct < 100; pct++) {
    std::vector<double> diff;
    for (auto& p : v)   {
      if (p.pct == pct) {
        double d = p.prs - p.expected;   // observed - expected
        if (p.pct >= 50) d = -d;         // flip sign for upper half
        diff.push_back(d);
        }

    }
    if (diff.size() < 2) continue;

    //double m = mean(diff);
    //bins[pct].effect = m;   // Python effect = mean residual
    //bins[pct].pval = t_pvalue(t, diff.size() - 1); 
    
    double m = mean(diff);
    double s = sd_pop(diff);  // ddof=0
    double se = s / std::sqrt((double)diff.size());
    double t = m / se;



                            //
    bins[pct].effect = m;  // Python effect = mean residual
    bins[pct].pval   = t_pvalue(t, (int)diff.size() - 1); 
    bins[pct].se     = se;
    bins[pct].eLo = m - 1.96 * se;
    bins[pct].eHi = m + 1.96 * se;
    
    //double m = mean(diff);
    //double s = sd(diff);
    //double t = m / (s / std::sqrt(diff.size()));
    //bins[pct].effect = m / s;
    //bins[pct].pval = t_pvalue(t, diff.size() - 1);
  }
}





//void PopData::calculate_qc() {
//        qc_ok = qc_middle_bins(bins, qc_min_p, qc_min_bin, qc_bonf_p);
//}







void PopData::calculate_bin_results(int regRange, int statRange) {
  // 1) run regression on chosen window; fills p.expected for ALL points
  run_reg(regRange);

  // 2) reset bins (always keep size 100)
  bins.assign(100, BinStats());

  // 3) choose which bins to compute stats for
  int stat_lo = (100 - statRange) / 2;
  int stat_hi = 99 - stat_lo;   // inclusive

  // 4) compute means for bins in stat range
  for (auto& p : v) {
    if (p.pct < stat_lo || p.pct > stat_hi) continue;
    auto& bn = bins[p.pct];
    bn.n++;
    bn.mean_prs += p.prs;
    bn.mean_exp += p.expected;
  }

  for (int pct = stat_lo; pct <= stat_hi; ++pct) {
    auto& bn = bins[pct];
    if (bn.n > 0) {
      bn.mean_prs /= bn.n;
      bn.mean_exp /= bn.n;
    }
  }

  // 5) compute per-bin p-values and SAVE the minimum
  min_pval_ = 1.0;
  min_pbin_ = -1;

  for (int pct = stat_lo; pct <= stat_hi; ++pct) {
    std::vector<double> diff;
    diff.reserve((size_t)bins[pct].n);

    for (auto& p : v) {
      if (p.pct == pct) {
        double d = p.prs - p.expected;   // keep YOUR existing residual definition
        if (p.pct >= 50) d = -d;         // keep YOUR existing sign flip
        diff.push_back(d);
      }
    }

    if (diff.size() < 2) continue;

    double m  = mean(diff);
    double s  = sd_pop(diff);                 // keep YOUR current choice
    double se = s / std::sqrt((double)diff.size());
    double t  = m / se;

    bins[pct].effect = m;
    bins[pct].pval   = t_pvalue(t, (int)diff.size() - 1);
    bins[pct].se     = se;
    bins[pct].eLo    = m - 1.96 * se;
    bins[pct].eHi    = m + 1.96 * se;

    if (bins[pct].pval < min_pval_) {
      min_pval_ = bins[pct].pval;
      min_pbin_ = pct;
    }
  }
}







void PopData::print_popout(int lo, int hi) const {
  for (int p : {lo, hi}) {
    const auto& b = bins[p];
    std::cout << "name" << name_ << " " << p
              << " n=" << b.n
              << " effect=" << b.effect
              << " p=" << b.pval
              << "\n";
  }
}






















void PopData::print_tail_header(std::ostream& out) {
  out
    << std::left  << std::setw(35) << "---"
    << std::right
    << std::setw(10) << "tail"
    << std::setw(10) << "len"
    << std::setw(10) << "r2"
    << std::setw(12) << "pv1"
    << std::setw(12) << "pv2"
    << std::setw(12) << "se1"
    << std::setw(12) << "se2"
    << std::setw(10) << "e1"
    << std::setw(10) << "e1Lo"
    << std::setw(10) << "e1Hi"
    << std::setw(10) << "e2"
    << std::setw(10) << "e2Lo"
    << std::setw(10) << "e2Hi"
    << std::setw(8)  << "QC"
    << "\n";
}




void PopData::print_tail_popout(std::ostream& out, bool qc_ok) const {
  const int lo = 0;
  const int hi = static_cast<int>(bins.size()) - 1;

  const auto& b1 = bins[lo];
  const auto& b2 = bins[hi];

  //const int len = b1.n + b2.n;
  const int len = static_cast<int>(v.size());
  const double r2 = 0.0;   // placeholder
  const bool qc = qc_ok;    

  //std::cerr << "QC min-p = " << qc_min_p << "\n"; 
  //std::cerr << "QC min-p loc = " << qc_min_bin << "\n"; 
  //std::cerr << "QC Bonf p = " << qc_bonf_p << "\n"; 
  //std::cerr << "QC min-p = " << min_pval() << "\n";
  //std::cerr << "QC min-p loc = " << min_pbin() << "\n";
  //std::cerr << "QC Bonf p = " << (min_pval() * 80.0) << "\n";

  out
    << std::left  << std::setw(35) << name_
    << std::right
    << std::setw(10) << tailSize_
    //<< std::setw(10) << std::fixed << std::setprecision(1) << 1.0
    << std::setw(10) << std::fixed << std::setprecision(2) << len
    << std::setw(10) << std::fixed << std::setprecision(4) << r2_ 
    << std::setw(12) << std::scientific << std::setprecision(4) << b1.pval
    << std::setw(12) << std::scientific << std::setprecision(4) << b2.pval
    << std::setw(12) << std::fixed << std::setprecision(3) << b1.se
    << std::setw(12) << std::fixed << std::setprecision(3) << b2.se
    << std::setw(10) << std::fixed << std::setprecision(4) << b1.effect
    << std::setw(10) << std::fixed << std::setprecision(3) << b1.eLo
    << std::setw(10) << std::fixed << std::setprecision(3) << b1.eHi
    << std::setw(10) << std::fixed << std::setprecision(4) << b2.effect
    << std::setw(10) << std::fixed << std::setprecision(3) << b2.eLo
    << std::setw(10) << std::fixed << std::setprecision(3) << b2.eHi
    << std::setw(8)  << std::boolalpha << qc
    << "\n";
}




