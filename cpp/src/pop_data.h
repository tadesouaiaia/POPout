#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <ostream>
#include <cmath>
#include <limits>


struct PopPoint {
  double prs;
  double pheno;
  int rank;
  int pct;
  double expected;
};

struct BinStats {
  int n = 0;
  double mean_prs = 0;
  double mean_exp = 0;
  double effect = 0;
  double pval = 1;
  double se  = NAN;
  double eLo = NAN;
  double eHi = NAN;
};

class PopData {
public:
    
    
    
  PopData(const std::string& path, const std::string& name, float tailSize = 1.0f, int regRange = 100);
  // NEW overloads (do not remove the old ones)
  // (optional) accessors if you want
  float tailSize() const { return tailSize_; }
  int   regRange()  const { return regRange_;  }


  void run_reg(int width);        // e.g. 100 = full, 90 = 5–95
  void calculate_means();

  void calculate_bin_results(int regRange = 100, int statRange = 100);
  double min_pval() const { return min_pval_; }
  int    min_pbin() const { return min_pbin_; }


  
  void print_popout(int lo, int hi) const;
  //static void print_tail_header(); 
  //void print_tail_popout() const;
  static void print_tail_header(std::ostream& out);
  void print_tail_popout(std::ostream& out, bool qc_ok = true) const;


  void calculate_qc();
  bool qc_ok = true;   // or declare + set in constructor if you prefer
  //double qc_bonf_p = 1.0;
  //double qc_min_p = 1.0;
  //int qc_min_bin = -1;











private:
  std::string name_;
  std::string path_;
  float tailSize_; 
  int regRange_;
  double r2_ = std::numeric_limits<double>::quiet_NaN(); 
  std::vector<PopPoint> v;
  std::vector<BinStats> bins;
  double min_pval_ = 1.0;
  int    min_pbin_ = -1;
  double a = 0; // intercept
  double b = 0; // slope
  void sort_and_rank();
};
