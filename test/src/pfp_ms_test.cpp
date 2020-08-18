/* pfp - prefix free parsing sa data structure test
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
   \file sa_test.cpp
   \brief sa_test.cpp build and test prefix-free parsing sa data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/


#include<iostream>
#include<vector>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>

#include <common.hpp>
#include <gtest/gtest.h>
#include <ms_w.hpp>
#include <pfp_ms_w.hpp>
#include <sdsl_ms_w.hpp>

extern "C" {
    #include<gsacak.h>
}

// Snippet from https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework
namespace testing
{
namespace internal
{
enum GTestColor
{
    COLOR_DEFAULT,
    COLOR_RED,
    COLOR_GREEN,
    COLOR_YELLOW
};

extern void ColoredPrintf(GTestColor color, const char *fmt, ...);
} // namespace internal
} // namespace testing
#define PRINTF(...)                                                                        \
    do                                                                                     \
    {                                                                                      \
        testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); \
        testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__);    \
    } while (0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
    ~TestCout()
    {
        PRINTF("%s", str().c_str());
    }
};

#define TEST_COUT TestCout()

//*************************************************************************************

std::string test_file;

typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<>>>, sdsl::lcp_support_sada<>> sdsl_cst_t;
typedef sdsl_cst_t::node_type sdsl_node_t;

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

    static const Names All[] = {L1, L2, L3};

    std::map<Names, std::pair<size_t, std::string>> Lengths = {
        {L1, {10, "l-10"}},
        {L2, {100, "l-100"}},
        {L3, {1000, "l-1000"}}};
} // namespace Query



class PFP_CST_Test : public ::testing::Test
{
protected:

    // Per-test-suite set-up.
    // Called before the first test in this test suite.
    // Can be omitted if not needed.
    static void SetUpTestCase()
    {
        TEST_COUT << "SetUpTestCase" << std::endl;
        // Load PFP_CST
        std::cout << "Loading Thresholds" << std::endl;

        ms = new ms_pointers<>(test_file);

        std::cout << "Building random access support" << std::endl;
        std::cout << "ALERT!!! w is hardcoded to be 10!" << std::endl;
        size_t w = 10;
        ra = new pfp_ra(test_file, w);

        // Load SDSL
        std::cout << "Loading SDSL CST" << std::endl;
        std::string filename = test_file + ".sdsl.cst";
        sdsl_cst = new sdsl_cst_t();
        sdsl::load_from_file(*sdsl_cst, filename);

        // Get wrappers
        pfp_w = new pfp_ms_w<ms_pointers<>, pfp_ra>(ms, ra);
        sdsl_w = new sdsl_ms_w<sdsl_cst_t>(sdsl_cst);

        // Generate samples
        TEST_COUT << "Generating samples" << std::endl;
        samples = new samples_t();
        generate_samples(1000, 0, samples, ra);
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    // Can be omitted if not needed.
    static void TearDownTestCase()
    {
        TEST_COUT << "TearDownTestCase" << std::endl;
        delete ra;
        delete ms;
        delete sdsl_cst;
        delete samples;
        delete pfp_w;
        delete sdsl_w;

        ra = nullptr;
        ms = nullptr;
        sdsl_cst = nullptr;
        samples = nullptr;
        pfp_w = nullptr;
        sdsl_w = nullptr;
    }

    // You can define per-test set-up logic as usual.
    virtual void SetUp() 
    {
        // TEST_COUT << "SetUp" << std::endl;
    }

    // You can define per-test tear-down logic as usual.
    virtual void TearDown() {}

    static void generate_samples(size_t Max_Sampling_Size, size_t seed, samples_t *samples, pfp_ra *ra)
    {
        srand(seed);
        size_t &n = ra->n;

        for (auto q : Query::All)
        {
            auto &l = Query::Lengths[q].first;
            (*samples).push_back(std::vector<query_t>());

            for (size_t i = 0; i < Max_Sampling_Size; ++i)
            {

                (*samples)[q].push_back(query_t());

                size_t pos = (rand() % n);
                size_t len = (rand() % l);

                while ((*samples)[q][i].size() < l)
                {
                    if (len == 0 || pos >= n)
                    {
                        pos = (rand() % n);
                        len = (rand() % l);
                    }

                    uint8_t c = ra->charAt(pos);
                    (*samples)[q][i].push_back(c);
                    len--;
                    pos++;
                }
            }
        }
    }

    // Some expensive resource shared by all tests.
    static const size_t w = 10;
    static ms_pointers<>* ms;
    static pfp_ra* ra;
    static sdsl_cst_t* sdsl_cst;
    static samples_t* samples;
    static ms_w *pfp_w;
    static ms_w *sdsl_w;
};

ms_pointers<> *PFP_CST_Test::ms = nullptr;
pfp_ra *PFP_CST_Test::ra = nullptr;
sdsl_cst_t *PFP_CST_Test::sdsl_cst = nullptr;
samples_t *PFP_CST_Test::samples = nullptr;
ms_w *PFP_CST_Test::pfp_w = nullptr;
ms_w *PFP_CST_Test::sdsl_w = nullptr;




TEST_F(PFP_CST_Test, CHILD)
{
    for (auto q : Query::All) // Move this as a parameter
    {
        size_t j = 0;
        for (const auto &query : (*samples)[q])
        {
            auto sdsl_res = sdsl_w->matching_statistics(query);
            auto pfp_res = pfp_w->matching_statistics(query);

            EXPECT_EQ(pfp_res.first.size(), sdsl_res.first.size());
            EXPECT_EQ(pfp_res.second.size(), sdsl_res.second.size());
            EXPECT_EQ(pfp_res.first.size(), sdsl_res.second.size());
            size_t n = pfp_res.first.size();
            for (size_t i = 0; i < n; ++i)
            {
                // EXPECT_EQ(pfp_res.first[i], sdsl_res.first[i]) << "At position: " << i;
                EXPECT_EQ(pfp_res.second[i], sdsl_res.second[i]) << "At position: " << i << " At position: " << j;
            }
            ++j;
        }
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " test_file " << std::endl;
        std::cout << " (1) Generates the SA, ISA and LCP;" << std::endl;
        std::cout << " (2) Generates LCE data structure and checks the result." << std::endl;
        return 1;
    }
    test_file = argv[1];
   
    return RUN_ALL_TESTS();
}

