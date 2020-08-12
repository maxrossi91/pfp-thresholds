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
   \file ms_w.hpp
   \brief ms_w.hpp ms wrapper crtp base.
   \author Massimiliano Rossi
   \date 11/08/2020
*/


#ifndef _MS_W_HH
#define _MS_W_HH

#include <iostream>
#include <sdsl/suffix_trees.hpp>

#include <benchmark/benchmark.h>

class ms_w{
public :
  typedef std::pair<std::vector<size_t>, std::vector<size_t>> ms_t;

  // The matchin statistics
  virtual ms_t matching_statistics(const std::vector<uint8_t>& pattern) = 0;

};

#endif /* end of include guard: _MS_W_HH */
