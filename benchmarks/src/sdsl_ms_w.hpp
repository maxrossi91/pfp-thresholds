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
   \file sdsl_ms_w.hpp
   \brief sdsl_ms_w.hpp cst wrapper crtp for sdsl.
   \author Massimiliano Rossi
   \date 08/08/2020
*/


#ifndef _SDSL_MS_W_HH
#define _SDSL_MS_W_HH

#include <iostream>
#include <ms_w.hpp>

#include <sdsl/io.hpp>
#include <sdsl/suffix_trees.hpp>

#include <benchmark/benchmark.h>

template<typename cst_t>
class sdsl_ms_w : public ms_w{
public :
  typedef typename ms_w::ms_t ms_t;

  cst_t* cst;


  sdsl_ms_w(cst_t* cst_):cst(cst_)
  {
    //Ntd
  }

  // The matchin statistics
  ms_t matching_statistics(const std::vector<uint8_t> &pattern)
  {
    typedef typename cst_t::size_type size_type;
    typedef typename cst_t::node_type node_type;

    size_t n2 = pattern.size();

    std::vector<size_t> pointers(n2);
    std::vector<size_t> lengths(n2);

    size_type cnt = 0;
    // sdsl::write_R_output("cst", "mstats", "begin", n2, cnt);
    size_type q = 0;                         // current match length
    size_type p2 = n2 - 1;                   // position in S2
    size_type i = 0, j = cst->csa.size() - 1; // \f$ \epsilon \f$ matches all suffixes of S1
    while (p2 + 1 > 0)
    {
      size_type lb, rb;
      // perform backward search on interval \f$ [i,j] \f$
      size_type size = sdsl::backward_search(cst->csa, i, j, pattern[p2], lb, rb);
      if (size > 0)
      {
        q = q + 1;
        i = lb;
        j = rb;
        p2 = p2 - 1;

        pointers[p2 + 1] = cst->sn(cst->node(i, i));
        lengths[p2 + 1] = q;
      }
      else if (i == 0 and j == cst->csa.size())
      {
        p2 = p2 - 1;
      }
      else
      {
        // map interval to a node of the cst and calculate parent
        node_type p = cst->parent(cst->node(i, j));
        q = cst->depth(p); // update match length
        i = cst->lb(p);    // update left bound
        j = cst->rb(p);    // update right bound
      }
      cnt += q;
    }

    return {pointers, lengths};
  }

};

#endif /* end of include guard: _SDSL_MS_W_HH */
