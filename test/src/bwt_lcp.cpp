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

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
    // space = thresholds.size() * sizeof(thresholds[0]);
    // space += thresholds_pos_s.size() * sizeof(thresholds_pos_s[0]);
    verbose("LCP size (bytes): ", space);
  }

  if (args.store)
  {
    verbose("Storing the LCP to file");
    std::string outfile = args.filename + std::string(".bwt2lcp.lcp");
    FILE* lcp_file;
    if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
      error("open() file " + outfile + " failed");

    for (size_t i = 0; i < M8.LCP.size(); ++i)
      if (fwrite(&M8.LCP[i], THRBYTES, 1, lcp_file) != 1)
        error("LCP write error 1");

    fclose(lcp_file);

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
