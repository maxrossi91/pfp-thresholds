/* ms_rle_string - Extension of the r-index rle_string to compute matching statistics
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
   \file ms_rle_string.hpp
   \brief ms_rle_string.hpp Extension of the r-index rle_string to compute matching statistics.
   \author Massimiliano Rossi
   \date 10/07/2020
*/

#ifndef _MS_RLE_STRING_HH
#define _MS_RLE_STRING_HH

#include <common.hpp>

#include <rle_string.hpp>


template <
    class sparse_bitvector_t = ri::sparse_sd_vector, //predecessor structure storing run length
    class string_t = ri::huff_string                 //run heads
    >
class ms_rle_string : public ri::rle_string<sparse_bitvector_t, string_t>
{
    public:

    ms_rle_string() : 
        ri::rle_string<sparse_bitvector_t, string_t>()
    {
        //NtD
    }

    /*
     * constructor: build structure on the input string
     * \param input the input string without 0x0 bytes in it.
     * \param B block size. The main sparse bitvector has R/B bits set (R being number of runs)
     *
     */
    ms_rle_string(string &input, ulint B = 2) : 
        ri::rle_string<sparse_bitvector_t, string_t>(input, B)
    {
        // NtD
    }

    ms_rle_string(std::ifstream &ifs, ulint B = 2) : 
        ri::rle_string<sparse_bitvector_t, string_t>(ifs, B)
    {

    }

    size_t number_of_runs_of_letter(uint8_t c)
    {
        return this->runs_per_letter[c].number_of_1();
    }

    size_t number_of_letter(uint8_t c)
    {
        return this->runs_per_letter[c].size();
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    ulint serialize(std::ostream &out) 
    {
        return ri::rle_string<sparse_bitvector_t, string_t>::serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {
        ri::rle_string<sparse_bitvector_t, string_t>::load(in);
    }

    private :
};

typedef ms_rle_string<ri::sparse_sd_vector> ms_rle_string_sd;
typedef ms_rle_string<ri::sparse_hyb_vector> ms_rle_string_hyb;

#endif /* end of include guard: _MS_RLE_STRING_HH */
