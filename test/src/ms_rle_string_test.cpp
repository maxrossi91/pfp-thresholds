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

class MS_RLE_String_Test : public ::testing::Test
{
protected:

    // Per-test-suite set-up.
    // Called before the first test in this test suite.
    // Can be omitted if not needed.
    static void SetUpTestCase()
    {
        TEST_COUT << "SetUpTestCase" << std::endl;
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    // Can be omitted if not needed.
    static void TearDownTestCase()
    {
        TEST_COUT << "TearDownTestCase" << std::endl;
    }

    // You can define per-test set-up logic as usual.
    virtual void SetUp() 
    {
        // TEST_COUT << "SetUp" << std::endl;
    }

    // You can define per-test tear-down logic as usual.
    virtual void TearDown() {}
};



TEST_F(MS_RLE_String_Test, RLBWT)
{
    std::string bwt_fname = test_file + ".bwt";

    TEST_COUT << "Construction from plain bwt" << std::endl;
    std::ifstream ifs(bwt_fname);
    ms_rle_string<> bwt(ifs);

    std::string bwt_heads_fname = bwt_fname + ".heads";
    std::ifstream ifs_heads(bwt_heads_fname);
    std::string bwt_len_fname = bwt_fname + ".len";
    std::ifstream ifs_len(bwt_len_fname);

    TEST_COUT << "Construction from rle bwt" << std::endl;
    ms_rle_string<> bwt1(ifs_heads, ifs_len);

    EXPECT_EQ(bwt.size(), bwt1.size());
    EXPECT_EQ(bwt.number_of_runs(), bwt1.number_of_runs());

    std::string bwt_string = bwt.toString();
    std::string bwt1_string = bwt1.toString();

    for (size_t i = 0; i < bwt_string.size(); ++i)
    {
        EXPECT_EQ(bwt_string[i], bwt1_string[i]) << "At position: " << i;
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

