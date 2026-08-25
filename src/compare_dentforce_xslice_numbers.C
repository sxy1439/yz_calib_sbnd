#include <TFile.h>
#include <TH2.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

const int nplanes = 3;
const int nxsections = 20;

const double xslice_low = -200.0;
const double xslice_width = 20.0;


// struct RatioSummary {
//   TString hist_name;
//   int plane;
//   int tpc;
//   int xsec;
//   double xlow;
//   double xhigh;

//   double global_median_original = 0.0;
//   double global_median_dentforce = 0.0;
//   double global_median_ratio = 0.0;

//   double median_bin_ratio = 0.0;
//   double mean_bin_ratio = 0.0;
//   double rms_bin_ratio = 0.0;

//   double min_bin_ratio = 0.0;
//   double max_bin_ratio = 0.0;

//   int n_bins_used = 0;
//   int n_original_nonzero = 0;
//   int n_dentforce_nonzero = 0;
//   int n_common = 0;
//   int n_original_only = 0;

//   double fraction_retained = 0.0;
//   double fraction_removed = 0.0;

//   double frac_outside_5pct = 0.0;
//   double frac_outside_10pct = 0.0;
// };

struct RatioSummary {
  TString hist_name;
  int plane;
  int tpc;
  int xsec;
  double xlow;
  double xhigh;

  // Medians using all nonzero bins in each histogram
  double median_original_all = 0.0;
  double median_dentforce_all = 0.0;
  double median_all_ratio = 0.0;

  // Medians using only common nonzero bins
  double median_original_common = 0.0;
  double median_dentforce_common = 0.0;
  double median_common_ratio = 0.0;

  // Bin-by-bin ratios, only where both original and DENTFORCE are nonzero
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
  int n_dentforce_only = 0;

  double fraction_retained = 0.0;
  double fraction_removed = 0.0;

  double frac_outside_5pct = 0.0;
  double frac_outside_10pct = 0.0;
};

double get_median(std::vector<double> vals)
{
  if (vals.empty()) return 0.0;

  std::sort(vals.begin(), vals.end());

  int n = vals.size();

  if (n % 2 == 1) {
    return vals[n / 2];
  }

  return 0.5 * (vals[n / 2 - 1] + vals[n / 2]);
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

// RatioSummary compare_one_hist(TH2* h_orig,
//                               TH2* h_dent,
//                               const TString& hname,
//                               int plane,
//                               int tpc,
//                               int xsec,
//                               double xlow,
//                               double xhigh)
// {
//   RatioSummary s;

//   s.hist_name = hname;
//   s.plane = plane;
//   s.tpc = tpc;
//   s.xsec = xsec;
//   s.xlow = xlow;
//   s.xhigh = xhigh;

//   if (!h_orig || !h_dent) {
//     std::cerr << "Missing histogram for " << hname << std::endl;
//     return s;
//   }

//   if (h_orig->GetNbinsX() != h_dent->GetNbinsX() ||
//       h_orig->GetNbinsY() != h_dent->GetNbinsY()) {
//     std::cerr << "Binning mismatch for " << hname << std::endl;
//     return s;
//   }

//   std::vector<double> orig_vals;
//   std::vector<double> dent_vals;
//   std::vector<double> ratios;

//   int n_outside_5 = 0;
//   int n_outside_10 = 0;

//   for (int ix = 1; ix <= h_orig->GetNbinsX(); ix++) {
//     for (int iy = 1; iy <= h_orig->GetNbinsY(); iy++) {

//       double orig = h_orig->GetBinContent(ix, iy);
//       double dent = h_dent->GetBinContent(ix, iy);

//       bool orig_valid = orig > 0.0;
//       bool dent_valid = dent > 0.0;

//       if (orig_valid) s.n_original_nonzero++;
//       if (dent_valid) s.n_dentforce_nonzero++;

//       if (orig_valid && dent_valid) {
//         s.n_common++;

//         orig_vals.push_back(orig);
//         dent_vals.push_back(dent);

//         double r = dent / orig;
//         ratios.push_back(r);

//         if (std::abs(r - 1.0) > 0.05) n_outside_5++;
//         if (std::abs(r - 1.0) > 0.10) n_outside_10++;
//       }

//       if (orig_valid && !dent_valid) {
//         s.n_original_only++;
//       }
//     }
//   }

//   s.n_bins_used = ratios.size();

//   if (ratios.empty()) {
//     return s;
//   }

//   s.global_median_original = get_median(orig_vals);
//   s.global_median_dentforce = get_median(dent_vals);

//   if (s.global_median_original > 0.0) {
//     s.global_median_ratio =
//       s.global_median_dentforce / s.global_median_original;
//   }

//   s.median_bin_ratio = get_median(ratios);

//   double sum = 0.0;
//   double sum2 = 0.0;

//   s.min_bin_ratio = ratios[0];
//   s.max_bin_ratio = ratios[0];

//   for (double r : ratios) {
//     sum += r;
//     sum2 += r * r;

//     if (r < s.min_bin_ratio) s.min_bin_ratio = r;
//     if (r > s.max_bin_ratio) s.max_bin_ratio = r;
//   }

//   s.mean_bin_ratio = sum / ratios.size();

//   double variance =
//     sum2 / ratios.size() - s.mean_bin_ratio * s.mean_bin_ratio;

//   if (variance < 0.0) variance = 0.0;

//   s.rms_bin_ratio = std::sqrt(variance);

//   s.frac_outside_5pct =
//     double(n_outside_5) / ratios.size();

//   s.frac_outside_10pct =
//     double(n_outside_10) / ratios.size();

//   if (s.n_original_nonzero > 0) {
//     s.fraction_retained =
//       double(s.n_common) / double(s.n_original_nonzero);

//     s.fraction_removed =
//       double(s.n_original_only) / double(s.n_original_nonzero);
//   }

//   return s;
// }

RatioSummary compare_one_hist(TH2* h_orig,
                              TH2* h_dent,
                              const TString& hname,
                              int plane,
                              int tpc,
                              int xsec,
                              double xlow,
                              double xhigh)
{
  RatioSummary s;

  s.hist_name = hname;
  s.plane = plane;
  s.tpc = tpc;
  s.xsec = xsec;
  s.xlow = xlow;
  s.xhigh = xhigh;

  if (!h_orig || !h_dent) {
    std::cerr << "Missing histogram for " << hname << std::endl;
    return s;
  }

  if (h_orig->GetNbinsX() != h_dent->GetNbinsX() ||
      h_orig->GetNbinsY() != h_dent->GetNbinsY()) {
    std::cerr << "Binning mismatch for " << hname << std::endl;
    return s;
  }

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

      // All original nonzero bins
      if (orig_valid) {
        s.n_original_nonzero++;
        orig_all_vals.push_back(orig);
      }

      // All DENTFORCE nonzero bins
      if (dent_valid) {
        s.n_dentforce_nonzero++;
        dent_all_vals.push_back(dent);
      }

      // Common bins only
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

      if (!orig_valid && dent_valid) {
        s.n_dentforce_only++;
      }
    }
  }

  s.n_bins_used = ratios.size();

  // All-bin medians
  s.median_original_all = get_median(orig_all_vals);
  s.median_dentforce_all = get_median(dent_all_vals);

  if (s.median_original_all > 0.0) {
    s.median_all_ratio =
      s.median_dentforce_all / s.median_original_all;
  }

  // Common-bin medians
  s.median_original_common = get_median(orig_common_vals);
  s.median_dentforce_common = get_median(dent_common_vals);

  if (s.median_original_common > 0.0) {
    s.median_common_ratio =
      s.median_dentforce_common / s.median_original_common;
  }

  // Bin-by-bin ratio summaries
  if (!ratios.empty()) {
    s.median_bin_ratio = get_median(ratios);

    double sum = 0.0;
    double sum2 = 0.0;

    s.min_bin_ratio = ratios[0];
    s.max_bin_ratio = ratios[0];

    for (double r : ratios) {
      sum += r;
      sum2 += r * r;

      if (r < s.min_bin_ratio) s.min_bin_ratio = r;
      if (r > s.max_bin_ratio) s.max_bin_ratio = r;
    }

    s.mean_bin_ratio = sum / ratios.size();

    double variance =
      sum2 / ratios.size() - s.mean_bin_ratio * s.mean_bin_ratio;

    if (variance < 0.0) variance = 0.0;

    s.rms_bin_ratio = std::sqrt(variance);

    s.frac_outside_5pct =
      double(n_outside_5) / ratios.size();

    s.frac_outside_10pct =
      double(n_outside_10) / ratios.size();
  }

  if (s.n_original_nonzero > 0) {
    s.fraction_retained =
      double(s.n_common) / double(s.n_original_nonzero);

    s.fraction_removed =
      double(s.n_original_only) / double(s.n_original_nonzero);
  }

  return s;
}

void compare_dentforce_xslice_numbers(
    const char* original_file,
    const char* dentforce_file,
    const char* out_csv = "../plots/dentforce_xslice_ratio_summary.csv",
    const char* hist_prefix = "CzyHist")
{
  TFile* f_orig = TFile::Open(original_file, "READ");
  TFile* f_dent = TFile::Open(dentforce_file, "READ");

  if (!f_orig || f_orig->IsZombie()) {
    std::cerr << "Could not open original file: "
              << original_file << std::endl;
    return;
  }

  if (!f_dent || f_dent->IsZombie()) {
    std::cerr << "Could not open DENTFORCE file: "
              << dentforce_file << std::endl;
    return;
  }

  std::ofstream csv(out_csv);

  // csv << "hist_type,"
  //     << "plane,"
  //     << "plane_label,"
  //     << "tpc,"
  //     << "tpc_label,"
  //     << "xsec,"
  //     << "xlow,"
  //     << "xhigh,"
  //     << "hist_name,"
  //     << "global_median_original,"
  //     << "global_median_dentforce,"
  //     << "global_median_ratio,"
  //     << "median_bin_ratio,"
  //     << "mean_bin_ratio,"
  //     << "rms_bin_ratio,"
  //     << "min_bin_ratio,"
  //     << "max_bin_ratio,"
  //     << "n_bins_used,"
  //     << "n_original_nonzero,"
  //     << "n_dentforce_nonzero,"
  //     << "n_common,"
  //     << "n_original_only,"
  //     << "fraction_retained,"
  //     << "fraction_removed,"
  //     << "frac_outside_5pct,"
  //     << "frac_outside_10pct"
  //     << "\n";

  csv << "hist_type,"
      << "plane,"
      << "plane_label,"
      << "tpc,"
      << "tpc_label,"
      << "xsec,"
      << "xlow,"
      << "xhigh,"
      << "hist_name,"
      << "median_original_all,"
      << "median_dentforce_all,"
      << "median_all_ratio,"
      << "median_original_common,"
      << "median_dentforce_common,"
      << "median_common_ratio,"
      << "median_bin_ratio,"
      << "mean_bin_ratio,"
      << "rms_bin_ratio,"
      << "min_bin_ratio,"
      << "max_bin_ratio,"
      << "n_bins_used,"
      << "n_original_nonzero,"
      << "n_dentforce_nonzero,"
      << "n_common,"
      << "n_original_only,"
      << "n_dentforce_only,"
      << "fraction_retained,"
      << "fraction_removed,"
      << "frac_outside_5pct,"
      << "frac_outside_10pct"
      << "\n";

  std::cout << "Comparing x-sliced "
            << hist_prefix
            << " histograms\n"
            << std::endl;

  for (int plane = 0; plane < nplanes; plane++) {
    for (int xs = 0; xs < nxsections; xs++) {

      int tpc = (xs < 10) ? 0 : 1;

      double xlow = xslice_low + xs * xslice_width;
      double xhigh = xlow + xslice_width;

      TString hname =
        Form("%s_%d_%d_xsec%d",
             hist_prefix,
             plane,
             tpc,
             xs);

      TH2* h_orig = (TH2*)f_orig->Get(hname);
      TH2* h_dent = (TH2*)f_dent->Get(hname);

      if (!h_orig || !h_dent) {
        std::cout << "Skipping missing histogram: "
                  << hname << std::endl;
        continue;
      }

      RatioSummary s =
        compare_one_hist(h_orig,
                         h_dent,
                         hname,
                         plane,
                         tpc,
                         xs,
                         xlow,
                         xhigh);

      std::cout << plane_label(plane)
                << " plane, "
                << tpc_label(tpc)
                << " TPC, "
                << "xsec" << xs
                << " (" << xlow << " < X < " << xhigh << " cm), "
                << hname << std::endl;

      // std::cout << "  global median original  = "
      //           << s.global_median_original << std::endl;

      // std::cout << "  global median DENTFORCE = "
      //           << s.global_median_dentforce << std::endl;

      // std::cout << "  global median ratio     = "
      //           << s.global_median_ratio << std::endl;

      // std::cout << "  median bin ratio        = "
      //           << s.median_bin_ratio << std::endl;

      // std::cout << "  mean bin ratio          = "
      //           << s.mean_bin_ratio << std::endl;

      // std::cout << "  RMS of bin ratio        = "
      //           << s.rms_bin_ratio << std::endl;

      // std::cout << "  bins used               = "
      //           << s.n_bins_used << std::endl;

      // std::cout << "  original nonzero bins   = "
      //           << s.n_original_nonzero << std::endl;

      // std::cout << "  DENTFORCE nonzero bins  = "
      //           << s.n_dentforce_nonzero << std::endl;

      // std::cout << "  fraction retained       = "
      //           << s.fraction_retained << std::endl;

      // std::cout << "  fraction removed        = "
      //           << s.fraction_removed << std::endl;

      // std::cout << "  fraction outside 5%     = "
      //           << s.frac_outside_5pct << std::endl;

      // std::cout << "  fraction outside 10%    = "
      //           << s.frac_outside_10pct << std::endl;

      // std::cout << std::endl;


      std::cout << "  median original all     = "
          << s.median_original_all << std::endl;

      std::cout << "  median DENTFORCE all    = "
                << s.median_dentforce_all << std::endl;

      std::cout << "  median all ratio        = "
                << s.median_all_ratio << std::endl;

      std::cout << "  median original common  = "
                << s.median_original_common << std::endl;

      std::cout << "  median DENTFORCE common = "
                << s.median_dentforce_common << std::endl;

      std::cout << "  median common ratio     = "
                << s.median_common_ratio << std::endl;

      std::cout << "  median bin ratio        = "
                << s.median_bin_ratio << std::endl;

      std::cout << "  mean bin ratio          = "
                << s.mean_bin_ratio << std::endl;

      std::cout << "  RMS of bin ratio        = "
                << s.rms_bin_ratio << std::endl;

      std::cout << "  min bin ratio           = "
                << s.min_bin_ratio << std::endl;

      std::cout << "  max bin ratio           = "
                << s.max_bin_ratio << std::endl;

      std::cout << "  bins used/common bins   = "
                << s.n_bins_used << std::endl;

      std::cout << "  original nonzero bins   = "
                << s.n_original_nonzero << std::endl;

      std::cout << "  DENTFORCE nonzero bins  = "
                << s.n_dentforce_nonzero << std::endl;

      std::cout << "  original-only bins      = "
                << s.n_original_only << std::endl;

      std::cout << "  DENTFORCE-only bins     = "
                << s.n_dentforce_only << std::endl;

      std::cout << "  fraction retained       = "
                << s.fraction_retained << std::endl;

      std::cout << "  fraction removed        = "
                << s.fraction_removed << std::endl;

      std::cout << "  fraction outside 5%     = "
                << s.frac_outside_5pct << std::endl;

      std::cout << "  fraction outside 10%    = "
                << s.frac_outside_10pct << std::endl;

      std::cout << std::endl;

      // csv << hist_prefix << ","
      //     << s.plane << ","
      //     << plane_label(s.plane) << ","
      //     << s.tpc << ","
      //     << tpc_label(s.tpc) << ","
      //     << s.xsec << ","
      //     << s.xlow << ","
      //     << s.xhigh << ","
      //     << s.hist_name << ","
      //     << s.global_median_original << ","
      //     << s.global_median_dentforce << ","
      //     << s.global_median_ratio << ","
      //     << s.median_bin_ratio << ","
      //     << s.mean_bin_ratio << ","
      //     << s.rms_bin_ratio << ","
      //     << s.min_bin_ratio << ","
      //     << s.max_bin_ratio << ","
      //     << s.n_bins_used << ","
      //     << s.n_original_nonzero << ","
      //     << s.n_dentforce_nonzero << ","
      //     << s.n_common << ","
      //     << s.n_original_only << ","
      //     << s.fraction_retained << ","
      //     << s.fraction_removed << ","
      //     << s.frac_outside_5pct << ","
      //     << s.frac_outside_10pct
      //     << "\n";


      csv << hist_prefix << ","
          << s.plane << ","
          << plane_label(s.plane) << ","
          << s.tpc << ","
          << tpc_label(s.tpc) << ","
          << s.xsec << ","
          << s.xlow << ","
          << s.xhigh << ","
          << s.hist_name << ","
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
          << s.n_dentforce_only << ","
          << s.fraction_retained << ","
          << s.fraction_removed << ","
          << s.frac_outside_5pct << ","
          << s.frac_outside_10pct
          << "\n";
    }
  }

  csv.close();

  f_orig->Close();
  f_dent->Close();

  std::cout << "Wrote CSV summary to: "
            << out_csv << std::endl;
}