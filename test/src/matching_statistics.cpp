/* matching_statistics - Computes the matching statistics from BWT and Thresholds
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
   \file matching_statistics.cpp
   \brief matching_statistics.cpp Computes the matching statistics from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 13/07/2020
*/

#define VERBOSE

#include <common.hpp>
#include <sdsl/io.hpp>
#include <ms_pointers.hpp>
#include <pfp_ra.hpp>
#include <malloc_count.h>

#include <iostream>
#include <omp.h>

struct Patterns
{
private:
  std::ifstream input_file;

public:
  using pattern_t = std::pair<std::string, std::vector<uint8_t>>;

  Patterns(std::string input_path) : input_file(input_path) {}

  pattern_t next();
  bool end();
};

Patterns::pattern_t Patterns::next()
{
  std::string header, read;
  std::vector<uint8_t> pattern;
  // header
  std::getline(input_file, header);
  // read
  std::getline(input_file, read);

  pattern_t out;
  out.first = header;
  std::copy(read.begin(), read.end(), std::back_inserter(out.second));
  return out;
}

bool Patterns::end()
{
  return input_file.eof();
}

int main(int argc, char *const argv[])
{

  Args args;
  parseArgs(argc, argv, args);

  // Building the r-index
  verbose("Building the matching statistics index");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  ms_pointers<> ms(args.filename);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Building random access");
  t_insert_start = std::chrono::high_resolution_clock::now();

  pfp_ra ra(args.filename, args.w);

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("Matching statistics index construction complete");
  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Memory peak: ", malloc_count_peak());
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(static)
  for (std::size_t i = 1; i <= 64; i++)
  {
    std::string file_path = args.patterns + "_" + std::to_string(i) + ".fa";

    std::ofstream f_pointers(file_path + ".pointers");
    std::ofstream f_lengths(file_path + ".lengths");

    verbose("Working on: ", file_path);

    if (!f_pointers.is_open())
      error("open() file " + std::string(file_path) + ".pointers failed");

    if (!f_lengths.is_open())
      error("open() file " + std::string(file_path) + ".lengths failed");

    Patterns patterns(file_path);

    while (not patterns.end())
    {
      Patterns::pattern_t pattern = patterns.next();

      auto pointers = ms.query(pattern.second);
      std::vector<size_t> lengths(pointers.size());
      size_t l = 0;

      for (size_t i = 0; i < pointers.size(); ++i)
      {
        size_t pos = pointers[i];
        while ((i + l) < pattern.second.size() && (pos + l) < ra.n && pattern.second[i + l] == ra.charAt(pos + l))
          ++l;

        lengths[i] = l;
        l = (l == 0 ? 0 : (l - 1));
      }

      f_pointers << pattern.first << endl;
      for (auto elem : pointers)
        f_pointers << elem << " ";
      f_pointers << endl;

      f_lengths << pattern.first << endl;
      for (auto elem : lengths)
        f_lengths << elem << " ";
      f_lengths << endl;
    }

    f_pointers.close();
    f_lengths.close();
  }

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
