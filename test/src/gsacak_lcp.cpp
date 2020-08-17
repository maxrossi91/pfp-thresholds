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
  if(args.is_fasta) 
    read_fasta_file(args.filename.c_str(), text);
  else
    read_file(args.filename.c_str(), text);
  verbose("Text size: ", text.size());

  // std::vector<char> tmp(args.w - 1, EndOfWord);
  // text.insert(text.begin(), tmp.begin(), tmp.end());

  uint8_t num_bytes = 1;

  verbose("Computing BWT and LCP");
  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

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

  
  lcp.erase(lcp.begin());

  // Opening output files
  FILE *lcp_file;
  std::string outfile = args.filename + std::string(".gsacak.lcp");
  if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");
    
  for(size_t i = 0; i < lcp.size(); ++i){
      size_t tmp_lcp = lcp[i];
      if (fwrite(&tmp_lcp, THRBYTES, 1, lcp_file) != 1)
        error("SA write error 1");
  }

  // Close output files
  fclose(lcp_file);

  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count();
  verbose("Elapsed time (s): ", time);

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
    verbose("Thresholds size (bytes): ", space);
  }

  if (args.store)
  {
    
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
  }
