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
   \file ptp_thresholds.cpp
   \brief ptp_thresholds.cpp test construction of MEM-thresholds from prefix-free parsing.
   \author Massimiliano Rossi
   \date 17/06/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_bwt.hpp>

#include <malloc_count.h>


int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);

  verbose("Reading text from file");
  std::vector<char> text;
  read_fasta_file(args.filename.c_str(), text);
  verbose("Text size: ", text.size());

  // std::vector<char> tmp(args.w - 1, EndOfWord);
  // text.insert(text.begin(), tmp.begin(), tmp.end());

  uint8_t num_bytes = 1;

  verbose("Computing BWT and LCP");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // Prepare Input
  sdsl::contains_no_zero_symbol(text, args.filename.c_str());
  sdsl::append_zero_symbol(text);
  // Construct SA
  sdsl::int_vector<> tmp_sa(text.size(), 0, sdsl::bits::hi(text.size()) + 1);
  sdsl::algorithm::calculate_sa((const unsigned char *)text.data(), text.size(), tmp_sa);
  std::vector<uint32_t> sa(tmp_sa.size());
  std::vector<uint32_t> isa(tmp_sa.size());
  std::vector<uint8_t> bwt(tmp_sa.size());
  for (size_t i = 0; i < tmp_sa.size(); ++i){
    sa[i] = tmp_sa[i];
    bwt[i] = text[(sa[i] == 0 ? sa.size() : sa[i]) - 1];
    isa[sa[i]] = i; 
  }

  std::vector<size_t> lcp(sa.size());

  LCP_array(&text[0], isa, sa, sa.size(), lcp);

  sdsl::rmq_succinct_sct<> rmq_lcp = sdsl::rmq_succinct_sct<>(&lcp);



  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("BWT and LCP construction complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  // auto time = std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();


  verbose("Building the thresholds");
  

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

  // This code gets timed

  // Opening output files
  FILE *thr_file;
  std::string outfile = args.filename + std::string(".sdsl.thr");
  if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  FILE *thr_pos_file;
  outfile = args.filename + std::string(".sdsl.thr_pos");
  if ((thr_pos_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  // std::vector<size_t> thresholds;
  // std::vector<size_t> thresholds_pos_s;
  std::vector<uint64_t> last_seen(256, 0);
  std::vector<bool> never_seen(256, true);


  never_seen[bwt[0]] = false;
  for(size_t i = 1; i < bwt.size(); ++i){
    if(bwt[i] == bwt[i-1])
      continue;
  
    if(never_seen[bwt[i]]){
      never_seen[bwt[i]] = false;
    }else{
      size_t j = rmq_lcp(last_seen[bwt[i]]+1, i);
      if (fwrite(&lcp[j], THRBYTES, 1, thr_file) != 1)
        error("SA write error 1");
      if (fwrite(&j, THRBYTES, 1, thr_pos_file) != 1)
        error("SA write error 1");
      // thresholds.push_back(lcp[j]);
      // thresholds_pos_s.push_back(j);
    }
    last_seen[bwt[i-1]] = i-1;
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
    // std::string outfile = args.filename + std::string(".sdsl.thr");
    // write_file(outfile.c_str(), thresholds);
    // outfile = args.filename + std::string(".sdsl.thr_pos");
    // write_file(outfile.c_str(), thresholds_pos_s);

    // verbose("Storing the BWT to file");
    // std::vector<uint8_t> my_bwt;
    // for(size_t i = 0; i < bwt.size(); ++i) my_bwt.push_back(bwt[i]);
    // outfile = args.filename + std::string(".sdsl.bwt");
    // write_file(outfile.c_str(), my_bwt);
    // verbose("Storing the LCP to file");
    // std::vector<int32_t> my_lcp;
    // for(size_t i = 0; i < lcp.size(); ++i ) my_lcp.push_back(lcp[i]);
    // outfile = args.filename + std::string(".sdsl.lcp");
    // write_file(outfile.c_str(), my_lcp);
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
  }
