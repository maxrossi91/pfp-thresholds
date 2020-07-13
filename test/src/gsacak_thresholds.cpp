/* gsaca_thresholds - Test of construction of MEM thresholds using gsacak
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
   \file gsacak_thresholds.cpp
   \brief gsacak_thresholds.cpp test construction of MEM-thresholds using gsacak.
   \author Massimiliano Rossi
   \date 26/06/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>

extern "C"
{
  #include <gsacak.h>
}

#include <malloc_count.h>

#define SABYTES 5


int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);

  verbose("Reading text from file");
  std::vector<uint8_t> text;
  read_fasta_file(args.filename.c_str(), text);
  verbose("Text size: ", text.size());

  // std::vector<char> tmp(args.w - 1, EndOfWord);
  // text.insert(text.begin(), tmp.begin(), tmp.end());

  uint8_t num_bytes = 1;

  verbose("Computing BWT and LCP");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // Prepare Input
  text.push_back(1);
  text.push_back(0);
  // Construct SA
  std::vector<uint_t> sa(text.size());
  std::vector<int_t> lcp(text.size());
  verbose("Computing SA, and LCP of the text");
  _elapsed_time(
    gsacak(&text[0], &sa[0], &lcp[0], nullptr, text.size())
  );

  verbose("Computing BWT of the text");
  // Opening output files
  FILE *ssa_file;
  std::string outfile = args.filename + std::string(".gsacak.ssa");
  if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  FILE *esa_file;
  outfile = args.filename + std::string(".gsacak.esa");
  if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  std::vector<uint8_t> bwt(text.size()-1);
  for (size_t i = 1; i < text.size(); ++i){
    bwt[i-1] = text[(sa[i] == 0 ? sa.size() : sa[i]) - 1];
    
    if(i > 1){
      if(bwt[i-1]!=bwt[i-2]){
        size_t ssa = sa[i];
        size_t esa = sa[i-1];
        // Output ssa and esa
        if (fwrite(&ssa, SABYTES, 1, ssa_file) != 1)
          error("SA write error 1");
        if (fwrite(&esa, SABYTES, 1, esa_file) != 1)
          error("SA write error 1");
      }
    }else{
      size_t ssa = sa[i];
      // Output ssa
      if (fwrite(&ssa, SABYTES, 1, ssa_file) != 1)
        error("SA write error 1");
    }
  }


  fclose(ssa_file);
  fclose(esa_file);

  outfile = args.filename + std::string(".gsacak.bwt");
  write_file(outfile.c_str(), bwt);

  lcp.erase(lcp.begin());

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
  outfile = args.filename + std::string(".gsacak.thr");
  if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  FILE *thr_pos_file;
  outfile = args.filename + std::string(".gsacak.thr_pos");
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
      // Write a zero so the positions of thresholds and BWT runs are the same
      size_t zero = 0;
      if (fwrite(&zero, THRBYTES, 1, thr_file) != 1)
        error("SA write error 1");
      if (fwrite(&zero, THRBYTES, 1, thr_pos_file) != 1)
        error("SA write error 1");
    }else{
      size_t j = rmq_lcp(last_seen[bwt[i]]+1, i);
      // size_t j = rmq_lcp(last_seen[bwt[i]]+1, i);
      size_t tmp_lcp = lcp[j];
      if (fwrite(&tmp_lcp, THRBYTES, 1, thr_file) != 1)
        error("SA write error 1");
      j--;
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
    // std::string outfile = args.filename + std::string(".gsacak.thr");
    // write_file(outfile.c_str(), thresholds);
    // outfile = args.filename + std::string(".gsacak.thr_pos");
    // write_file(outfile.c_str(), thresholds_pos_s);

    // verbose("Storing the BWT to file");
    // std::vector<uint8_t> my_bwt;
    // for(size_t i = 0; i < bwt.size(); ++i) my_bwt.push_back(bwt[i]);
    // outfile = args.filename + std::string(".gsacak.bwt");
    // write_file(outfile.c_str(), my_bwt);
    // verbose("Storing the LCP to file");
    // std::vector<int32_t> my_lcp;
    // for(size_t i = 0; i < lcp.size(); ++i ) my_lcp.push_back(lcp[i]);
    // outfile = args.filename + std::string(".gsacak.lcp");
    // write_file(outfile.c_str(), my_lcp);
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
  }
