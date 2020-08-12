/* pfp - prefix free parsing compress suffix tree
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
   \file pfp_ms_w.hpp
   \brief pfp_ms_w.hpp ms wrapper crtp for pfp_ms.
   \author Massimiliano Rossi
   \date 11/08/2020
*/


#ifndef _PFP_MS_W_HH
#define _PFP_MS_W_HH

#include <iostream>
#include <ms_w.hpp>

#include <ms_pointers.hpp>
#include <pfp_ra.hpp>

#include <benchmark/benchmark.h>

template<typename ms_p_t, typename ra_t>
class pfp_ms_w : public ms_w{
public :
  typedef ms_w::ms_t ms_t;
  
  ms_p_t *ms_p;
  ra_t *ra;

  pfp_ms_w(ms_p_t* ms_p_, ra_t* ra_):ms_p(ms_p_),ra(ra_)
  {
    //Ntd
  }

  // The matchin statistics
  ms_t matching_statistics(const std::vector<uint8_t> &pattern)
  {
    auto pointers = ms_p->query(pattern);
    std::vector<size_t> lengths(pointers.size());
    size_t l = 0;
    for (size_t i = 0; i < pointers.size(); ++i)
    {
      size_t pos = pointers[i];
      while ((i + l) < pattern.size() && (pos + l) < ra->n && pattern[i + l] == ra->charAt(pos + l))
        ++l;

      lengths[i] = l;
      l = (l == 0 ? 0 : (l - 1));
    }

    return {pointers,lengths};
  }

};

#endif /* end of include guard: _PFP_MS_W_HH */
