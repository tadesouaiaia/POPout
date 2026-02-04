#include "stats.h"

#include <cmath>
#include <limits>




// Regularized incomplete beta via continued fraction
static double betacf(double a, double b, double x) {
  const int MAXIT = 100;
  const double EPS = 3e-7;
  const double FPMIN = 1e-30;

  double qab = a + b;
  double qap = a + 1.0;
  double qam = a - 1.0;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  double h = d;

  for (int m = 1; m <= MAXIT; m++) {
    int m2 = 2 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (std::fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d * c;

    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (std::fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < EPS) break;
  }
  return h;
}

static double betai(double a, double b, double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;

  double bt = std::exp(
    std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
    + a * std::log(x) + b * std::log(1.0 - x)
  );

  if (x < (a + 1.0) / (a + b + 2.0))
    return bt * betacf(a, b, x) / a;
  else
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

// two-sided Student-t p-value
double t_pvalue(double t, int df) {
  double x = df / (df + t * t);
  double p = betai(0.5 * df, 0.5, x);
  return p;
}





void reg_with_intercept(const std::vector<double>& x,
                        const std::vector<double>& y,
                        double& a,
                        double& b) {
  int n = x.size();
  double sx=0, sy=0, sxx=0, sxy=0;

  for (int i = 0; i < n; i++) {
    sx += x[i];
    sy += y[i];
    sxx += x[i] * x[i];
    sxy += x[i] * y[i];
  }

  double mx = sx / n;
  double my = sy / n;

  b = (sxy - n * mx * my) / (sxx - n * mx * mx);
  a = my - b * mx;
}

double mean(const std::vector<double>& x) {
  double s = 0;
  for (double v : x) s += v;
  return s / x.size();
}

double sd(const std::vector<double>& x) {
  double m = mean(x);
  double s = 0;
  for (double v : x) s += (v - m) * (v - m);
  return std::sqrt(s / (x.size() - 1));
}

double sd_pop(const std::vector<double>& x) {
      double m = mean(x);
        double s = 0.0;
          for (double v : x) s += (v - m) * (v - m);
            return std::sqrt(s / (double)x.size());  // ddof=0
}                                   // }



// assumes mean(x) and sd_pop(x) already exist




double se_pop(const std::vector<double>& x) {
    const int n = (int)x.size();
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    return sd_pop(x) / std::sqrt((double)n);   // ddof=0, consistent with sd_pop
}

void ci95_mean_pop(const std::vector<double>& x, double& lo, double& hi) {
    lo = hi = std::numeric_limits<double>::quiet_NaN();
    const int n = (int)x.size();
    if (n < 2) return;

    const double m  = mean(x);
    const double se = se_pop(x);
    if (!std::isfinite(m) || !std::isfinite(se) || se == 0) return;

    constexpr double z = 1.96;   // normal 95% CI
    lo = m - z * se;
    hi = m + z * se;
}




double corr_sums(double n,
                 double sumx, double sumy,
                 double sumxx, double sumyy,
                 double sumxy)
{
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();

    const double sxx = n * sumxx - sumx * sumx;
    const double syy = n * sumyy - sumy * sumy;
    const double sxy = n * sumxy - sumx * sumy;

    if (sxx <= 0 || syy <= 0)
        return std::numeric_limits<double>::quiet_NaN();

    return sxy / std::sqrt(sxx * syy);
}





bool qc_middle_bins(const std::vector<BinStats>& bins, double& min_p, int& min_bin, double& bonf_p) {
    int nbins = (int)bins.size();

    int start = 10;
    int end   = nbins - 10;   // exclusive
    int n_qc  = end - start;  // 80 for nbins=100

    min_p = 1.0;
    min_bin = -1;

    for (int pct = start; pct < end; ++pct) {
        double p = bins[pct].pval;
        if (p < min_p) {
            min_p = p;
            min_bin = pct;
        }
    }

    bonf_p = min_p * n_qc;
    return bonf_p >= 0.05;
}







