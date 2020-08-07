/* pfp_thresholds - Test of construction of MEM thresholds from prefix free parsing data structures
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
   \file pfp_thresholds.cpp
   \brief pfp_thresholds.cpp test construction of MEM-thresholds from prefix-free parsing.
   \author Massimiliano Rossi
   \date 17/06/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <pfp.hpp>
#include <pfp_lcp.hpp>

#include <malloc_count.h>

int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);

  // TODO: Include cmd line option to load from file/compute
  // verbose("Loading PFP data structures from file");
  // std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // pf_parsing<> pf;
  // std::string filename = args.filename + pf.filesuffix();
  // sdsl::load_from_file(pf, filename);

  // Computing prefix-free parsing
  verbose("Window size set to: " , args.w);

  verbose("Computing PFP data structures");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  pf_parsing pf(args.filename, args.w);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PFP DS construction complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  // auto time = std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();

  verbose("Memory peak: ", malloc_count_peak());
  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

  // This code gets timed

  verbose("Building the LCP");

  
  pfp_lcp lcp(pf);

  
  // verbose("Memory peak: ", malloc_count_peak());
  // verbose("Building the thresholds - constructing thresholds");

  // // Opening output files
  // FILE *thr_file;
  // std::string outfile = args.filename + std::string(".thr");
  // if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
  //   error("open() file " + outfile + " failed");
  
  // FILE *thr_pos_file;
  // outfile = args.filename + std::string(".thr_pos");
  // if ((thr_pos_file = fopen(outfile.c_str(), "w")) == nullptr)
  //   error("open() file " + outfile + " failed");

  // // std::vector<size_t> thresholds;
  // // std::vector<size_t> thresholds_pos_s;
  // std::vector<uint64_t> last_seen(256, 0);
  // std::vector<bool> never_seen(256, true);

  // sdsl::rmq_succinct_sct<> rmq_min_s = sdsl::rmq_succinct_sct<>(&lcp.min_s);

  // for(size_t i = 1; i < lcp.heads.size(); ++i){
  //   if (never_seen[lcp.heads[i]])
  //   {
  //     never_seen[lcp.heads[i]] = false;

  //     // Write a zero so the positions of thresholds and BWT runs are the same
  //     size_t zero = 0;
  //     if (fwrite(&zero, THRBYTES, 1, thr_file) != 1)
  //       error("SA write error 1");
  //     if (fwrite(&zero, THRBYTES, 1, thr_pos_file) != 1)
  //       error("SA write error 1");
  //   }
  //   else
  //   {
  //     size_t j = rmq_min_s(last_seen[lcp.heads[i]] + 1, i - 1);
  //     if (fwrite(&lcp.min_s[j], THRBYTES, 1, thr_file) != 1)
  //       error("SA write error 1");
  //     if (fwrite(&lcp.pos_s[j], THRBYTES, 1, thr_pos_file) != 1)
  //       error("SA write error 1");
  //     // thresholds.push_back(min_s[j]);
  //     // thresholds_pos_s.push_back(pos_s[j]);
  //   }
  //   last_seen[lcp.heads[i]] = i;
  // }

  // // Close output files
  // fclose(thr_file);
  // fclose(thr_pos_file);

  // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  // auto time = std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count();
  // verbose("Elapsed time (s): ", time);

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
    lcp.pos_s.erase(lcp.pos_s.begin());
    lcp.min_s.erase(lcp.min_s.begin());
    verbose("Storing the LCP to file");
    std::string outfile = args.filename + std::string(".pfp.lcp.pos");
    write_file(outfile.c_str(), lcp.pos_s);
    outfile = args.filename + std::string(".pfp.lcp.min");
    write_file(outfile.c_str(), lcp.min_s);
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
  }
