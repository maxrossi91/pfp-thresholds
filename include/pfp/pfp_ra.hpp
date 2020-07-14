/* pfp - prefix free parsing for random access
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
   \file pfp.hpp
   \brief pfp.hpp define and build the prefix-free parsing data structures for random access.
   \author Massimiliano Rossi
   \date 13/07/2020
   \note This is a short version of the prefix free parsing in https://github.com/maxrossi91/pfp-data-structures
*/

#ifndef _PFP_RA_HH
#define _PFP_RA_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>


class pfp_ra{
public:
  // Dictionary
  std::vector<uint8_t> d;

  sdsl::bit_vector b_d; // Starting position of each phrase in D
  sdsl::bit_vector::rank_1_type rank_b_d;
  sdsl::bit_vector::select_1_type select_b_d;

  // Parse
  std::vector<uint32_t> p;

  size_t n; // Size of the text
  size_t w; // Size of the window

  sdsl::bit_vector b_p;
  sdsl::bit_vector::rank_1_type rank_b_p;
  sdsl::bit_vector::select_1_type select_b_p;

  typedef size_t size_type;

  // Default constructor for load
  pfp_ra() {}

  pfp_ra(std::vector<uint8_t> &d_,
             std::vector<uint32_t> &p_,
             size_t w_) : 
            d(d_),
            p(p_),
            w(w_)
  {


    verbose("Building b_d");
    _elapsed_time(build_b_d());

    // Compute the length of the string;
    compute_n();

    verbose("Building b_p");
    _elapsed_time(build_b_p());

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  pfp_ra( std::string filename, size_t w_):
              w(w_)
  {
    read_dictionary(filename);
    read_parse(filename);

    verbose("Building b_d");
    _elapsed_time(build_b_d());

    // Compute the length of the string;
    compute_n();

    verbose("Building b_p");
    _elapsed_time(build_b_p());

    // print_sizes();

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  inline size_t length_of_phrase(size_t id)
  {
    assert(id > 0);
    return select_b_d(id + 1) - select_b_d(id) - 1; // to remove the EndOfWord
  }

  inline size_t n_phrases()
  {
    return rank_b_d(d.size() - 1);
  }

  void print_sizes()
  {

    verbose("Parse");
    verbose("Size of p: ", p.size() * sizeof(p[0]));

    verbose("Dictionary");
    verbose("Size of d: ", d.size() * sizeof(d[0]));
    verbose("Size of b_d: ", sdsl::size_in_bytes(b_d));
    verbose("Size of rank_b_d: ", sdsl::size_in_bytes(rank_b_d));
    verbose("Size of select_b_d: ", sdsl::size_in_bytes(select_b_d));

    verbose("PFP");
    verbose("Size of b_p: ", sdsl::size_in_bytes(b_p));
    verbose("Size of rank_b_p: ", sdsl::size_in_bytes(rank_b_p));
    verbose("Size of select_b_p: ", sdsl::size_in_bytes(select_b_p));

  }

  
  void clear_unnecessary_elements(){
    // Reducing memory tentative

    
  }

  uint8_t charAt(
      uint64_t pos //!< 0-based position
  ) const
  {
    // This is because, the text is considered to be cyclic.
    pos = (pos + w) % n;

    assert(pos < n);

    size_t phrase_idx = rank_b_p(pos + 1) - 1;
    uint32_t phrase_id = p[phrase_idx];
    assert(pos >= select_b_p(phrase_idx + 1));
    size_t offset = pos - select_b_p(phrase_idx + 1);

    size_t pos_in_d = select_b_d(phrase_id) + offset;
    return d[pos_in_d];
  }

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += my_serialize(p, out, child, "p");
    written_bytes += my_serialize(d, out, child, "d");
    written_bytes += b_d.serialize(out, child, "b_d");
    written_bytes += rank_b_d.serialize(out, child, "rank_b_d");
    written_bytes += select_b_d.serialize(out, child, "select_b_d");
    written_bytes += b_p.serialize(out, child, "b_p");
    written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    written_bytes += select_b_p.serialize(out, child, "select_b_p");
    written_bytes += sdsl::write_member(n, out, child, "n");
    written_bytes += sdsl::write_member(w, out, child, "w");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    my_load(p, in);
    my_load(d, in);
    b_d.load(in);
    rank_b_d.load(in);
    select_b_d.load(in);
    b_p.load(in);
    rank_b_p.load(in);
    select_b_p.load(in);
    sdsl::read_member(n, in);
    sdsl::read_member(w, in);
  }

  std::string filesuffix() const
  {
    return ".pf.ds.ra";
  }

  protected:
    void read_dictionary(std::string filename)
    {
      // Building dictionary from file
      std::string tmp_filename = filename + std::string(".dict");
      read_file(tmp_filename.c_str(), d);
      assert(d[0] == Dollar);

      // Prepending w dollars to d
      // 1. Count how many dollars there are
      int i = 0;
      int n_dollars = 0;
      while (i < d.size() && d[i++] == Dollar)
        ++n_dollars;
      std::vector<uint8_t> dollars(w - n_dollars, Dollar);
      d.insert(d.begin(), dollars.begin(), dollars.end());
    }

    void read_parse(std::string filename)
    {
      // Building dictionary from file
      std::string tmp_filename = filename + std::string(".parse");
      read_file(tmp_filename.c_str(), p);
      p.push_back(0); // this is the terminator for the sacak algorithm
    }

    void build_b_d()
    {
      // Building the bitvector with a 1 in each starting position of each phrase in D
      b_d.resize(d.size());
      for (size_t i = 0; i < b_d.size(); ++i)
        b_d[i] = false; // bug in resize
      b_d[0] = true;    // Mark the first phrase
      for (size_t i = 1; i < d.size(); ++i)
        b_d[i] = (d[i - 1] == EndOfWord);
      b_d[d.size() - 1] = true; // This is necessary to get the length of the last phrase

      rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
      select_b_d = sdsl::bit_vector::select_1_type(&b_d);
    }

    void build_b_p()
    {
      // Build the bitvector storing the position of the beginning of each phrase.
      b_p.resize(this->n); // all should be initialized at false by sdsl
      for (size_t i = 0; i < b_p.size(); ++i)
        b_p[i] = false; // bug in resize
      b_p[0] = true;    // phrase_0 becomes phrase 1

      size_t i = 0;

      for (int j = 0; j < p.size() - 2; ++j)
      { // -2 because the beginning of the last phrase is in position 0
        // p[i]: phrase_id
        assert(p[j] != 0);
        // phrase_length: select_b_d(p[i]+1)-select_b_d(p[i]);
        i += length_of_phrase(p[j]) - w;
        b_p[i] = true;
      }

      // Build rank and select on Sp
      rank_b_p = sdsl::bit_vector::rank_1_type(&b_p);
      select_b_p = sdsl::bit_vector::select_1_type(&b_p);
    }

    void compute_n()
    {
      // Compute the length of the string;
      n = 0;
      for (int j = 0; j < p.size() - 1; ++j)
      {
        // p[j]: phrase_id
        assert(p[j] != 0);
        n += length_of_phrase(p[j]) - w;
      }
    }
};

#endif /* end of include guard: _PFP_RA_HH */
