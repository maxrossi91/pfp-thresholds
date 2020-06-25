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

template <typename T>
void write_file(const char *filename, std::vector<T> &ptr)
{
    struct stat filestat;
    FILE *fd;

    if ((fd = fopen(filename, "w")) == nullptr)
        error("open() file " + std::string(filename) + " failed");

    size_t length = ptr.size();
    if ((fwrite(&ptr[0], sizeof(T), length, fd)) != length)
        error("fwrite() file " + std::string(filename) + " failed");

    fclose(fd);
}

int main(int argc, char *const argv[])
{

    Args args;
    parseArgs(argc, argv, args);

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
    verbose("Reading text from file");
    std::vector<char> text;
    read_fasta_file(args.filename.c_str(), text);
    verbose("Text size: ", text.size());

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count();
    verbose("Elapsed time (s): ", time);

    auto mem_peak = malloc_count_peak();
    verbose("Memory peak: ", malloc_count_peak());

    size_t space = 0;
    if (args.memo)
    {
        space = text.size() * sizeof(text[0]);
        verbose("Text size (bytes): ", space);
    }

    if (args.store)
    {
        verbose("Storing Text to file");
        std::string outfile = args.filename + std::string(".plain");
        write_file(outfile.c_str(), text);
    }

    if (args.csv)
        std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

    return 0;
}
