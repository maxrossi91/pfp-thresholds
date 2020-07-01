/* pfp - prefix free parsing 
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
   \brief pfp.hpp define and build the prefix-free parsing data structures.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note This is a short version of the prefix free parsing in https://github.com/maxrossi91/pfp-data-structures
*/

#ifndef _PFP_HH
#define _PFP_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}

#include<dictionary.hpp>
#include<parse.hpp>

class pf_parsing{
public:
  
  dictionary dict;
  parse pars;
  std::vector<uint32_t> freq;
  size_t n; // Size of the text
  size_t w; // Size of the window

  sdsl::bit_vector b_p;
  sdsl::bit_vector::rank_1_type rank_b_p;
  sdsl::bit_vector::select_1_type select_b_p;

  typedef size_t size_type;

  // Default constructor for load
  pf_parsing() {}

  pf_parsing(std::vector<uint8_t> &d_,
             std::vector<uint32_t> &p_,
             std::vector<uint32_t> &freq_,
             size_t w_) : 
            dict(d_, w_),
            pars(p_, dict.n_phrases() + 1),
            freq(freq_),
            w(w_)
  {
    // Uploading the frequency file
    assert(freq[0] == 0);

    // Compute the length of the string;
    compute_n();

    verbose("Computing b_p");
    _elapsed_time(compute_b_p());

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  pf_parsing( std::string filename, size_t w_):
              dict(filename, w_),
              pars(filename,dict.n_phrases()+1),
              freq(1,0),
              w(w_)
  {
    // Uploading the frequency file
    uint32_t *occ;
    size_t d_words;
    std::string tmp_filename = filename + std::string(".occ");
    read_file<uint32_t> (tmp_filename.c_str(), occ, d_words);
    freq.insert(freq.end(),occ, occ + d_words);


    // Compute the length of the string;
    compute_n();

    // b_p(pfp.n,0);
    verbose("Computing b_p");
    _elapsed_time(compute_b_p());

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  void compute_b_p() {
    // Build the bitvector storing the position of the beginning of each phrase.
    b_p.resize(this->n); // all should be initialized at false by sdsl
    for(size_t i = 0; i < b_p.size(); ++i)
      b_p[i] = false; // bug in resize
    b_p[0] = true; // phrase_0 becomes phrase 1
    
    size_t i = 0;
    
    for(int j = 0; j < pars.p.size()-2; ++j){ // -2 because the beginning of the last phrase is in position 0
      // p[i]: phrase_id
      assert(pars.p[j] != 0);
      // phrase_length: select_b_d(p[i]+1)-select_b_d(p[i]);
      i += dict.length_of_phrase(pars.p[j]) - w;
      b_p[i] = true;
    }

    // Build rank and select on Sp
    rank_b_p = sdsl::bit_vector::rank_1_type(&b_p);
    select_b_p = sdsl::bit_vector::select_1_type(&b_p);
  }

  void compute_n(){
    // Compute the length of the string;
    n = 0;
    for (int j = 0; j < pars.p.size() - 1; ++j)
    {
      // parse.p[j]: phrase_id
      assert(pars.p[j] != 0);
      n += dict.length_of_phrase(pars.p[j]) - w;
    }
    //n += w; // + w because n is the length including the last w markers
    //n += w - 1; // Changed after changind b_d in dict // -1 is for the first dollar + w because n is the length including the last w markers
  }


  void clear_unnecessary_elements(){
    // NtD
    pars.isaP.clear(); pars.isaP.shrink_to_fit();
  }

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += dict.serialize(out, child, "dictionary");
    written_bytes += pars.serialize(out, child, "parse");
    written_bytes += my_serialize(freq, out, child, "frequencies");
    written_bytes += sdsl::write_member(n, out, child, "n");
    written_bytes += sdsl::write_member(w, out, child, "w");
    written_bytes += b_p.serialize(out, child, "b_p");
    written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    written_bytes += select_b_p.serialize(out, child, "select_b_p");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    dict.load(in);
    pars.load(in);
    my_load(freq, in);
    sdsl::read_member(n, in);
    sdsl::read_member(w, in);
    b_p.load(in);
    rank_b_p.load(in, &b_p);
    select_b_p.load(in, &b_p);
  }

  std::string filesuffix() const
  {
    return ".pf.ds.thr";
  }


};

#endif /* end of include guard: _PFP_HH */
