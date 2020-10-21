/* sdsl_matching_statistics - Computes the matching statistics using sdsl
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
   \file sdsl_matching_statistics.cpp
   \brief sdsl_matching_statistics.cpp Computes the matching statistics using sdsl.
   \author Massimiliano Rossi
   \date 14/07/2020
   \note from http://simongog.github.io/sdsl/a00266_source.html
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>


#include <sdsl/io.hpp>
#include <sdsl/suffix_trees.hpp>


#include <malloc_count.h>

// ************************* From: http://simongog.github.io/sdsl/a00266_source.html
template <class Cst>
std::vector<size_t> cst_matching_statistics(const Cst &cst, unsigned char *S2, typename Cst::size_type n2)
{
  typedef typename Cst::size_type size_type;
  typedef typename Cst::node_type node_type;

  std::vector<size_t> lengths(n2);

  size_type cnt = 0;
  // sdsl::write_R_output("cst", "mstats", "begin", n2, cnt);
  size_type q = 0;                         // current match length
  size_type p2 = n2 - 1;                   // position in S2
  size_type i = 0, j = cst.csa.size() - 1; // \f$ \epsilon \f$ matches all suffixes of S1
  while (p2 + 1 > 0)
  {
    size_type lb, rb;
    // perform backward search on interval \f$ [i,j] \f$
    size_type size = sdsl::backward_search(cst.csa, i, j, S2[p2], lb, rb);
    if (size > 0)
    {
      q = q + 1;
      i = lb;
      j = rb;
      p2 = p2 - 1;

      lengths[p2] = q;
    }
    else if (i == 0 and j == cst.csa.size())
    {
      p2 = p2 - 1;
    }
    else
    {
      // map interval to a node of the cst and calculate parent
      node_type p = cst.parent(cst.node(i, j));
      q = cst.depth(p); // update match length
      i = cst.lb(p);    // update left bound
      j = cst.rb(p);    // update right bound
    }
    cnt += q;
  }
  // sdsl::write_R_output("cst", "mstats", "end", n2, cnt);
  return lengths;
}

//********************************************************************************

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;

std::vector<pattern_t> read_patterns(std::string filename)
{
  // Open File
  FILE *fd;
  if ((fd = fopen(filename.c_str(), "r")) == nullptr)
    error("open() file " + filename + " failed");

  std::vector<pattern_t> patterns;

  pattern_t pattern;

  char c;
  while (fread(&c, sizeof(char), 1, fd) == 1)
  {
    if (c == '>')
    {
      if(pattern.second.size()>0)
        patterns.push_back(pattern);

      pattern.first.clear();
      pattern.second.clear();
      
      pattern.first.append(1,c);
      while (fread(&c, sizeof(char), 1, fd) == 1 && c != '\n')
        pattern.first.append(1,c);
    }
    else
    {
      pattern.second.push_back(c);
      while (fread(&c, sizeof(char), 1, fd) == 1 && c != '\n')
        pattern.second.push_back(c);
    }
  }

  if (pattern.second.size() > 0)
    patterns.push_back(pattern);
    
  fclose(fd);

  return patterns;
}

int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index

  verbose("Building the matching statistics index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  verbose("Reading text from file");
  std::vector<char> text;
  read_fasta_file(args.filename.c_str(), text);
  verbose("Text size: ", text.size());

  uint8_t num_bytes = 1;

  verbose("Computing CST of the text");
  sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> cst;
  auto time = _elapsed_time(
      sdsl::construct_im(cst, static_cast<const char *>(&text[0]), num_bytes);
  );

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  
  verbose("Reading patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  std::vector<pattern_t> patterns = read_patterns(args.patterns);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  
  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

  // std::ofstream f_pointers(args.patterns + ".sdsl.pointers");
  std::ofstream f_lengths(args.patterns + ".sdsl.lengths");

  // if (!f_pointers.is_open())
  //     error("open() file " + std::string(args.filename) + ".sdsl.pointers failed");

  if (!f_lengths.is_open())
      error("open() file " + std::string(args.filename) + ".sdsl.lengths failed");

  for(auto pattern: patterns)
  {

    std::vector<size_t> lengths = cst_matching_statistics(cst, &pattern.second[0], pattern.second.size());
    // auto pointers = ms.query(pattern.second);
    // std::vector<size_t> lengths(pointers.size());
    // size_t l = 0;
    // for(size_t i = 0; i < pointers.size(); ++i)
    // {
    //   size_t pos = pointers[i];
    //   while ((i + l) < pattern.second.size() && (pos + l) < ra.n && pattern.second[i + l] == ra.charAt(pos + l))
    //     ++l;

    //   lengths[i] = l;
    //   l = (l == 0? 0: (l-1));
    // }


    // f_pointers << pattern.first << endl;
    // for(auto elem: pointers)
    //   f_pointers << elem << " ";
    // f_pointers << endl;

    f_lengths << pattern.first << std::endl;
    for(auto elem: lengths)
      f_lengths << elem << " ";
    f_lengths << std::endl;
    
  }
  
  // f_pointers.close();
  f_lengths.close();

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());






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
