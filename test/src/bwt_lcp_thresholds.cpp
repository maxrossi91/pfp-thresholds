/* bwt_lcp_thresholds - Test of construction of MEM thresholds from bwt and rlbwt2lcp
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file bwt_lcp_thresholds.cpp
   \brief bwt_lcp_thresholds.cpp test construction of MEM-thresholds from bwt and rlbwt2lcp.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note Adapted main from https://github.com/nicolaprezza/rlbwt2lcp
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>

#include <lcp.hpp>
#include <unistd.h>
#include <dna_bwt_n.hpp>

#include <malloc_count.h>


template<typename bwt_t>
void build_thresholds(Args args){

  char TERM = EndOfDict;

  verbose("Alphabet: A,C,G,T,'" , TERM , "'");
  verbose("Loading and indexing BWT ... ");

  bwt_t bwt = bwt_t(args.filename, TERM);

  verbose("Done. Size of BWT: " , bwt.size());

  verbose("LCP construction");

  std::chrono::high_resolution_clock::time_point t_lcp_start = std::chrono::high_resolution_clock::now();

  lcp<bwt_t, uint64_t> M8(&bwt);

  std::chrono::high_resolution_clock::time_point t_lcp_end = std::chrono::high_resolution_clock::now();
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_lcp_end - t_lcp_start).count());

  sdsl::rmq_succinct_sct<> rmq_lcp = sdsl::rmq_succinct_sct<>(&M8.LCP);

  verbose("Building the thresholds");

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

  // This code gets timed

  // Opening output files
  FILE *thr_file;
  std::string outfile = args.filename + std::string(".rlbwt2lcp.thr");
  if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  FILE *thr_pos_file;
  outfile = args.filename + std::string(".rlbwt2lcp.thr_pos");
  if ((thr_pos_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  // std::vector<size_t> thresholds;
  // std::vector<size_t> thresholds_pos_s;
  std::vector<size_t> last_seen(256, 0);
  std::vector<bool> never_seen(256, true);

  size_t r = 0;
  never_seen[bwt[0]] = false;
  for (size_t i = 1; i < bwt.size(); ++i)
  {
    if (bwt[i] == bwt[i - 1])
      continue;

    if (never_seen[bwt[i]])
    {
      never_seen[bwt[i]] = false;
    }
    else
    {
      size_t j = rmq_lcp(last_seen[bwt[i]] + 1, r);
      if (fwrite(&M8.LCP[j], THRBYTES, 1, thr_file) != 1)
        error("SA write error 1");
      if (fwrite(&M8.MIN[j], THRBYTES, 1, thr_pos_file) != 1)
        error("SA write error 1");
      // thresholds.push_back(M8.LCP[j]);
      // thresholds_pos_s.push_back(M8.MIN[j]);
    }
    last_seen[bwt[i - 1]] = r;
    r++;
  }

  // Close output files
  fclose(thr_file);
  fclose(thr_pos_file);

  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count();
  verbose("Elapsed time (s): ", time);

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
    // space = thresholds.size() * sizeof(thresholds[0]);
    // space += thresholds_pos_s.size() * sizeof(thresholds_pos_s[0]);
    verbose("Thresholds size (bytes): ", space);
  }

  if (args.store)
  {
    // verbose("Storing the Thresholds to file");
    // std::string outfile = args.filename + std::string(".rlbwt2lcp.thr");
    // write_file(outfile.c_str(), thresholds);
    // outfile = args.filename + std::string(".rlbwt2lcp.thr_pos");
    // write_file(outfile.c_str(), thresholds_pos_s);
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;
}

int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);

  bool containsN = false;


  containsN = hasN(args.filename);

  if (not containsN)
    build_thresholds<dna_bwt_t>(args);
  else
    build_thresholds<dna_bwt_n_t>(args);

  verbose("Done. ");

  return 0;
  }
