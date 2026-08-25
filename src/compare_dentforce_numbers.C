#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

double median(std::vector<double>& v)
{
  if (v.empty()) return 0.0;

  std::sort(v.begin(), v.end());

  int n = v.size();

  if (n % 2 == 1) {
    return v[n/2];
  }

  return 0.5 * (v[n/2 - 1] + v[n/2]);
}

// struct RatioSummary {
//   double global_median_original = 0.0;
//   double global_median_dentforce = 0.0;
//   double global_median_ratio = 0.0;

//   double median_bin_ratio = 0.0;
//   double mean_bin_ratio = 0.0;
//   double rms_bin_ratio = 0.0;
//   double min_bin_ratio = 0.0;
//   double max_bin_ratio = 0.0;

//   int n_bins_used = 0;
//   double frac_outside_5pct = 0.0;
//   double frac_outside_10pct = 0.0;
// };

// RatioSummary compare_hists(TH2* h_orig, TH2* h_dent)
// {
//   RatioSummary s;

//   std::vector<double> orig_vals;
//   std::vector<double> dent_vals;
//   std::vector<double> ratios;

//   int n_outside_5 = 0;
//   int n_outside_10 = 0;

//   for (int ix = 1; ix <= h_orig->GetNbinsX(); ix++) {
//     for (int iy = 1; iy <= h_orig->GetNbinsY(); iy++) {

//       double orig = h_orig->GetBinContent(ix, iy);
//       double dent = h_dent->GetBinContent(ix, iy);

//       // Just realised that this only counts bins that are non-zero in both histograms 
//       if (orig <= 0) continue;
//       if (dent <= 0) continue;

//       double r = dent / orig;

//       orig_vals.push_back(orig);
//       dent_vals.push_back(dent);
//       ratios.push_back(r);

//       if (std::abs(r - 1.0) > 0.05) n_outside_5++;
//       if (std::abs(r - 1.0) > 0.10) n_outside_10++;
//     }
//   }

//   s.n_bins_used = ratios.size();

//   if (ratios.empty()) return s;

//   s.global_median_original = median(orig_vals);
//   s.global_median_dentforce = median(dent_vals);

//   if (s.global_median_original != 0) {
//     s.global_median_ratio = s.global_median_dentforce / s.global_median_original;
//   }

//   s.median_bin_ratio = median(ratios);

//   double sum = 0.0;
//   double sumsq = 0.0;

//   s.min_bin_ratio = ratios[0];
//   s.max_bin_ratio = ratios[0];

//   for (double r : ratios) {
//     sum += r;
//     sumsq += r * r;

//     if (r < s.min_bin_ratio) s.min_bin_ratio = r;
//     if (r > s.max_bin_ratio) s.max_bin_ratio = r;
//   }

//   s.mean_bin_ratio = sum / ratios.size();

//   double mean_sq = sumsq / ratios.size();
//   s.rms_bin_ratio = std::sqrt(mean_sq - s.mean_bin_ratio * s.mean_bin_ratio);

//   s.frac_outside_5pct = double(n_outside_5) / ratios.size();
//   s.frac_outside_10pct = double(n_outside_10) / ratios.size();

//   return s;
// }


struct RatioSummary {
  double median_original_all = 0.0;
  double median_dentforce_all = 0.0;
  double median_all_ratio = 0.0;

  double median_original_common = 0.0;
  double median_dentforce_common = 0.0;
  double median_common_ratio = 0.0;

  double median_bin_ratio = 0.0;
  double mean_bin_ratio = 0.0;
  double rms_bin_ratio = 0.0;
  double min_bin_ratio = 0.0;
  double max_bin_ratio = 0.0;

  int n_bins_used = 0;
  int n_original_nonzero = 0;
  int n_dentforce_nonzero = 0;
  int n_common = 0;
  int n_original_only = 0;

  double fraction_retained = 0.0;
  double fraction_removed = 0.0;

  double frac_outside_5pct = 0.0;
  double frac_outside_10pct = 0.0;
};

RatioSummary compare_hists(TH2* h_orig, TH2* h_dent)
{
  RatioSummary s;

  std::vector<double> orig_all_vals;
  std::vector<double> dent_all_vals;

  std::vector<double> orig_common_vals;
  std::vector<double> dent_common_vals;

  std::vector<double> ratios;

  int n_outside_5 = 0;
  int n_outside_10 = 0;

  for (int ix = 1; ix <= h_orig->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= h_orig->GetNbinsY(); iy++) {

      double orig = h_orig->GetBinContent(ix, iy);
      double dent = h_dent->GetBinContent(ix, iy);

      bool orig_valid = orig > 0.0;
      bool dent_valid = dent > 0.0;

      if (orig_valid) {
        s.n_original_nonzero++;
        orig_all_vals.push_back(orig);
      }

      if (dent_valid) {
        s.n_dentforce_nonzero++;
        dent_all_vals.push_back(dent);
      }

      if (orig_valid && dent_valid) {
        s.n_common++;

        orig_common_vals.push_back(orig);
        dent_common_vals.push_back(dent);

        double r = dent / orig;
        ratios.push_back(r);

        if (std::abs(r - 1.0) > 0.05) n_outside_5++;
        if (std::abs(r - 1.0) > 0.10) n_outside_10++;
      }

      if (orig_valid && !dent_valid) {
        s.n_original_only++;
      }
    }
  }

  s.n_bins_used = ratios.size();

  s.median_original_all = median(orig_all_vals);
  s.median_dentforce_all = median(dent_all_vals);

  if (s.median_original_all != 0.0) {
    s.median_all_ratio = s.median_dentforce_all / s.median_original_all;
  }

  s.median_original_common = median(orig_common_vals);
  s.median_dentforce_common = median(dent_common_vals);

  if (s.median_original_common != 0.0) {
    s.median_common_ratio = s.median_dentforce_common / s.median_original_common;
  }

  if (ratios.empty()) return s;

  s.median_bin_ratio = median(ratios);

  double sum = 0.0;
  double sumsq = 0.0;

  s.min_bin_ratio = ratios[0];
  s.max_bin_ratio = ratios[0];

  for (double r : ratios) {
    sum += r;
    sumsq += r * r;

    if (r < s.min_bin_ratio) s.min_bin_ratio = r;
    if (r > s.max_bin_ratio) s.max_bin_ratio = r;
  }

  s.mean_bin_ratio = sum / ratios.size();

  double mean_sq = sumsq / ratios.size();
  double variance = mean_sq - s.mean_bin_ratio * s.mean_bin_ratio;
  if (variance < 0.0) variance = 0.0;

  s.rms_bin_ratio = std::sqrt(variance);

  s.frac_outside_5pct = double(n_outside_5) / ratios.size();
  s.frac_outside_10pct = double(n_outside_10) / ratios.size();

  if (s.n_original_nonzero > 0) {
    s.fraction_retained = double(s.n_common) / s.n_original_nonzero;
    s.fraction_removed = double(s.n_original_only) / s.n_original_nonzero;
  }

  return s;
}
const char* plane_label(int plane)
{
  if (plane == 0) return "U";
  if (plane == 1) return "V";
  if (plane == 2) return "Y";
  return "Unknown";
}

const char* tpc_label(int tpc)
{
  if (tpc == 0) return "East";
  if (tpc == 1) return "West";
  return "Unknown";
}

void compare_dentforce_numbers(const char* original_file,
                               const char* dentforce_file,
                               const char* out_csv = "../plots/dentforce_ratio_summary.csv",
                               const char* hist_prefix = "zy")
{
  TFile* f_orig = TFile::Open(original_file, "READ");
  TFile* f_dent = TFile::Open(dentforce_file, "READ");

  if (!f_orig || f_orig->IsZombie()) {
    std::cerr << "ERROR: Could not open original file: " << original_file << std::endl;
    return;
  }

  if (!f_dent || f_dent->IsZombie()) {
    std::cerr << "ERROR: Could not open DENTFORCE file: " << dentforce_file << std::endl;
    return;
  }

  std::ofstream csv(out_csv);


  csv << "hist_type,plane,tpc,hist_name,"
    << "median_original_all,median_dentforce_all,median_all_ratio,"
    << "median_original_common,median_dentforce_common,median_common_ratio,"
    << "median_bin_ratio,mean_bin_ratio,rms_bin_ratio,min_bin_ratio,max_bin_ratio,"
    << "n_bins_used,n_original_nonzero,n_dentforce_nonzero,n_common,n_original_only,"
    << "fraction_retained,fraction_removed,"
    << "frac_outside_5pct,frac_outside_10pct\n";

  std::cout << "\nComparing " << hist_prefix << " histograms\n\n";

  for (int plane = 0; plane < 3; plane++) {
    for (int tpc = 0; tpc < 2; tpc++) {

      TString hname = Form("%s_%i_%i", hist_prefix, plane, tpc);

      TH2* h_orig = (TH2*)f_orig->Get(hname);
      TH2* h_dent = (TH2*)f_dent->Get(hname);

      if (!h_orig || !h_dent) {
        std::cout << "Skipping missing hist: " << hname << std::endl;
        continue;
      }

      RatioSummary s = compare_hists(h_orig, h_dent);

      std::cout << plane_label(plane) << " plane, "
                << tpc_label(tpc) << " TPC, "
                << hname << "\n";

      // std::cout << "  global median original  = " << s.global_median_original << "\n";
      // std::cout << "  global median DENTFORCE = " << s.global_median_dentforce << "\n";
      // std::cout << "  global median ratio     = " << s.global_median_ratio << "\n";
      // std::cout << "  median bin ratio        = " << s.median_bin_ratio << "\n";
      // std::cout << "  mean bin ratio          = " << s.mean_bin_ratio << "\n";
      // std::cout << "  RMS of bin ratio        = " << s.rms_bin_ratio << "\n";
      // std::cout << "  bins used               = " << s.n_bins_used << "\n";
      // std::cout << "  fraction outside 5%     = " << s.frac_outside_5pct << "\n";
      // std::cout << "  fraction outside 10%    = " << s.frac_outside_10pct << "\n\n";




      std::cout << "  median original all     = " << s.median_original_all << "\n";
      std::cout << "  median DENTFORCE all    = " << s.median_dentforce_all << "\n";
      std::cout << "  median all ratio        = " << s.median_all_ratio << "\n";

      std::cout << "  median original common  = " << s.median_original_common << "\n";
      std::cout << "  median DENTFORCE common = " << s.median_dentforce_common << "\n";
      std::cout << "  median common ratio     = " << s.median_common_ratio << "\n";

      std::cout << "  median bin ratio        = " << s.median_bin_ratio << "\n";
      std::cout << "  mean bin ratio          = " << s.mean_bin_ratio << "\n";
      std::cout << "  RMS of bin ratio        = " << s.rms_bin_ratio << "\n";

      std::cout << "  min bin ratio           = " << s.min_bin_ratio << "\n";
      std::cout << "  max bin ratio           = " << s.max_bin_ratio << "\n";

      std::cout << "  bins used/common bins   = " << s.n_bins_used << "\n";
      std::cout << "  original nonzero bins   = " << s.n_original_nonzero << "\n";
      std::cout << "  DENTFORCE nonzero bins  = " << s.n_dentforce_nonzero << "\n";
      std::cout << "  common bins             = " << s.n_common << "\n";
      std::cout << "  original-only bins      = " << s.n_original_only << "\n";

      std::cout << "  fraction retained       = " << s.fraction_retained << "\n";
      std::cout << "  fraction removed        = " << s.fraction_removed << "\n";

      std::cout << "  fraction outside 5%     = " << s.frac_outside_5pct << "\n";
      std::cout << "  fraction outside 10%    = " << s.frac_outside_10pct << "\n\n";

      csv << hist_prefix << ","
          << plane_label(plane) << ","
          << tpc_label(tpc) << ","
          << hname << ","
          << s.median_original_all << ","
          << s.median_dentforce_all << ","
          << s.median_all_ratio << ","
          << s.median_original_common << ","
          << s.median_dentforce_common << ","
          << s.median_common_ratio << ","
          << s.median_bin_ratio << ","
          << s.mean_bin_ratio << ","
          << s.rms_bin_ratio << ","
          << s.min_bin_ratio << ","
          << s.max_bin_ratio << ","
          << s.n_bins_used << ","
          << s.n_original_nonzero << ","
          << s.n_dentforce_nonzero << ","
          << s.n_common << ","
          << s.n_original_only << ","
          << s.fraction_retained << ","
          << s.fraction_removed << ","
          << s.frac_outside_5pct << ","
          << s.frac_outside_10pct << "\n";
    }
  }

  csv.close();

  f_orig->Close();
  f_dent->Close();

  std::cout << "Saved summary CSV to: " << out_csv << std::endl;
}