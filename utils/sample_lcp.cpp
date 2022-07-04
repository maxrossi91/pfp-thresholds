/* fasta2plain - Convert a fasta file to plain text.
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
   \file fasta2plain.cpp
   \brief fasta2plain.cpp convert a fasta file to plain text.
   \author Massimiliano Rossi
   \date 24/06/2020
*/

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <malloc_count.h>

int main(int argc, char *const argv[])
{

    Args args;
    parseArgs(argc, argv, args);

    // Building the r-index

    verbose("Reading the BWT and LCP array");
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    // Open output files
    std::string bwt_filename = args.filename + ".bwt";
    std::string lcp_filename = args.filename + ".lcp";
    std::string slcp_filename = args.filename + ".slcp";

    FILE *bwt;
    FILE *lcp;
    FILE *slcp;

    if ((bwt = fopen(bwt_filename.c_str(), "r")) == nullptr)
        error("open() file " + std::string(bwt_filename) + " failed");

    if ((lcp = fopen(lcp_filename.c_str(), "r")) == nullptr)
        error("open() file " + std::string(lcp_filename) + " failed");

    if ((slcp = fopen(slcp_filename.c_str(), "w")) == nullptr)
        error("open() file " + std::string(slcp_filename) + " failed");

    // Get the BWT length
    struct stat filestat;
    int fn = fileno(bwt);

    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(bwt_filename) + " failed");

    if (filestat.st_size % sizeof(uint8_t) != 0)
        error("invilid file " + std::string(bwt_filename));

    size_t n = filestat.st_size;

    // Start processing
    size_t curr_lcp;
    uint8_t prev, curr;
    size_t len = 1;

    if ((fread(&prev, sizeof(uint8_t), 1, bwt)) != 1)
        error("fread() file " + std::string(bwt_filename) + " failed");

    if ((fread((char*)&curr_lcp, SSABYTES, 1, lcp)) != 1)
        error("fread() file " + std::string(lcp_filename) + " failed");

    if ((fwrite((char*)&curr_lcp, SSABYTES, 1, slcp)) != 1)
        error("fwrite() file " + std::string(slcp_filename) + " failed");

    for (size_t i = 0; i < n - 1; ++i)
    {
        if ((fread(&curr, sizeof(uint8_t), 1, bwt)) != 1)
            error("fread() file " + std::string(bwt_filename) + " failed");

        if ((fread(&curr_lcp, SSABYTES, 1, lcp)) != 1)
            error("fread() file " + std::string(lcp_filename) + " failed");
        if (curr == prev)
        {
            len++;
        }
        else
        {
            if ((fwrite(&curr_lcp, SSABYTES, 1, slcp)) != 1)
                error("fwrite() file " + std::string(slcp_filename) + " failed");

            prev = curr;
            len = 1;
        }
    }


    fclose(bwt);
    fclose(lcp);
    fclose(slcp);

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

    auto mem_peak = malloc_count_peak();
    verbose("Memory peak: ", malloc_count_peak());

    size_t space = 0;
    if (args.memo)
    {
    }

    if (args.store)
    {
    }

    if (args.csv)
        std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

    return 0;
}
