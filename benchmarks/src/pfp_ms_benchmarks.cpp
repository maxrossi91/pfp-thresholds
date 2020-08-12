/* pfp - prefix free parsing lce structure test
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
   \file pfp_ms_benchmarks.cpp
   \brief pfp_ms_benchmarks.cpp load and test matching statistics.
   \author Massimiliano Rossi
   \date 11/08/2020
*/


#include<iostream>
#include<vector>

#include <sdsl/suffix_trees.hpp>

#include <pfp.hpp>
#include <ms_w.hpp>
#include <pfp_ms_w.hpp>
#include <sdsl_ms_w.hpp>


extern "C" {
    #include<gsacak.h>
}

#include <benchmark/benchmark.h>

typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst_t;

typedef std::vector<uint8_t> query_t;
typedef std::vector<std::vector<query_t>> samples_t;

namespace Query
{
    enum Names
    {
        L1,
        L2,
        L3
    };

    static const Names All[] = {L1,L2,L3};

    std::map<Names, std::pair<size_t, std::string>> Lengths = {
        {L1, {10, "l-10"}}, 
        {L2, {100, "l-100"}}, 
        {L3, {1000, "l-1000"}}};
} // namespace Query



void generate_samples(size_t Max_Sampling_Size, size_t seed, samples_t& samples, pfp_ra* ra)
{
    srand(seed);
    size_t& n = ra->n;

    for (auto q : Query::All)
    {
        auto& l = Query::Lengths[q].first;
        samples.push_back(std::vector<query_t>());

        for (size_t i = 0; i < Max_Sampling_Size; ++i)
        {

            samples[q].push_back(query_t());

            size_t pos = (rand() % n);
            size_t len = (rand() % l);

            while (samples[q][i].size() < l)
            {
                if(len == 0 || pos >= n)
                {
                    pos = (rand() % n);
                    len = (rand() % l);
                }

                uint8_t c = ra->charAt(pos);
                samples[q][i].push_back(c);
                len--; pos++;
            }
            
        }

    }

}

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State &_state)
{
    for (auto _ : _state)
    {
        std::string empty_string;
    }

    _state.counters["Size(bytes)"] = 0;
    _state.counters["Bits_x_Symbol"] = 0;
    _state.counters["Queries"] = 0;
    _state.counters["Time_x_Query"] = 0;
}
BENCHMARK(BM_WarmUp);

typedef std::function < void(benchmark::State &, ms_w *&, size_t &, samples_t &, size_t &, Query::Names& )> lambda_t;

// Benchmark Parent query
auto BM_MS =
[](benchmark::State &_state, auto _idx, const auto &_size, auto &_samples, const auto &_length, const auto& q ) {
    for (auto _ : _state)
    {
        for (const auto& query : _samples[q])
        {
            benchmark::DoNotOptimize(_idx->matching_statistics(query));
        }
    }

    _state.counters["Size(bytes)"] = _size;
    _state.counters["Bits_x_Symbol"] = _size * 8.0 / _length;
    _state.counters["Queries"] = _samples[q].size();
    _state.counters["Time_x_Query"] = benchmark::Counter(
        _samples[q].size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};


int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " test_file " << std::endl;
        return 1;
    }
    std::string test_file = argv[1];


    // std::vector<std::string> _data = {
    //     "../../../data/Chr19/chr19.1.fa",
    //     "../../../data/Chr19/chr19.2.fa",
    //     "../../../data/Chr19/chr19.4.fa",
    //     "../../../data/Chr19/chr19.8.fa",
    //     "../../../data/Chr19/chr19.16.fa",
    //     "../../../data/Chr19/chr19.32.fa",
    //     "../../../data/Chr19/chr19.64.fa",
    //     "../../../data/Chr19/chr19.128.fa",
    //     "../../../data/Chr19/chr19.256.fa",
    //     "../../../data/Chr19/chr19.512.fa",
    //     "../../../data/Chr19/chr19.1000.fa",
    //     "../../../data/Salmonella/salmonella.50.fa.fix",
    //     "../../../data/Salmonella/salmonella.100.fa.fix",
    //     "../../../data/Salmonella/salmonella.500.fa.fix",
    //     "../../../data/Salmonella/salmonella.1000.fa.fix",
    //     "../../../data/Salmonella/salmonella.5000.fa.fix",
    //     "../../../data/Salmonella/salmonella.10000.fa.fix",
    //     "../../../data/pizzachili/repcorpus/real/einstein.en.txt",
    //     "../../../data/pizzachili/repcorpus/real/world_leaders",
    //     "../../../data/pizzachili/repcorpus/real/cere"
    // }




    // Load PFP_CST
    std::cout << "Loading Thresholds"<< std::endl;

    ms_pointers<> ms(test_file);

    std::cout << "Building random access support"<< std::endl;
    std::cout << "ALERT!!! w is hardwired to be 10!"<< std::endl;
    size_t w = 10;
    pfp_ra ra(test_file, w);



    // Load SDSL
    std::cout << "Loading SDSL CST"<< std::endl;
    sdsl_cst_t sdsl_cst;
    std::string filename = test_file + ".sdsl.cst";
    sdsl::load_from_file(sdsl_cst, filename);

    // Saple setup
    size_t Max_Sampling_Size = 100;
    size_t seed = 0;
    size_t n = sdsl_cst.size();

    samples_t samples;

    std::cout << "Generating " << Max_Sampling_Size << " samples" << std::endl;
    generate_samples(Max_Sampling_Size, seed, samples, &ra);

    // Get wrappers
    ms_w* pfp_w = new pfp_ms_w<ms_pointers<>, pfp_ra>(&ms, &ra);
    ms_w* sdsl_w = new sdsl_ms_w<sdsl_cst_t>(&sdsl_cst);

    sdsl::nullstream ns;

    size_t pfp_size = ms.serialize(ns) + sdsl::size_in_bytes(ra);//sdsl::size_in_bytes(ms) + sdsl::size_in_bytes(ra);
    size_t sdsl_size = sdsl::size_in_bytes(sdsl_cst);

    std::cout << "Registering benchmarks" << std::endl;
    std::string pfp_s = "pfp";
    std::string sdsl_s = "sdsl";

    std::vector<std::pair<std::string, std::pair<ms_w*, size_t> > >  csts = {
        {"pfp", {pfp_w,pfp_size}},
        {"sdsl",{sdsl_w,sdsl_size}}};


    for(auto cst: csts){

        for (auto q : Query::All)
        {
            auto op = Query::Lengths[q];

            auto bm_name = cst.first + "-" + op.second;

            benchmark::RegisterBenchmark(bm_name.c_str(), BM_MS , cst.second.first, cst.second.second, samples, n, q);
        }
    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    for(auto cst: csts)
        delete cst.second.first;
    return 0;
}