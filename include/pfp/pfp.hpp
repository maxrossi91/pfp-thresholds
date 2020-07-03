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
  std::vector<uint_t> freq;
  size_t n; // Size of the text
  size_t w; // Size of the window

  std::vector<size_t> s_lcp_T; // LCP array of T sampled in corrispondence of the beginning of each phrase.
  sdsl::rmq_succinct_sct<> rmq_s_lcp_T;
  
  std::vector<int_t>  ilist;            // Inverted list of phrases of P in BWT_P
  sdsl::bit_vector ilist_s; // The ith 1 is in correspondence of the first occurrence of the ith phrase
  sdsl::bit_vector::select_1_type select_ilist_s;

  sdsl::bit_vector b_p;
  sdsl::bit_vector::rank_1_type rank_b_p;
  sdsl::bit_vector::select_1_type select_b_p;

  typedef size_t size_type;

  // Default constructor for load
  pf_parsing() {}

  pf_parsing(std::vector<uint8_t> &d_,
             std::vector<uint32_t> &p_,
             std::vector<uint_t> &freq_,
             size_t w_) : 
            dict(d_, w_),
            pars(p_, dict.n_phrases() + 1),
            freq(freq_),
            s_lcp_T(1,0),
            ilist(pars.p.size()),
            ilist_s(pars.p.size()+ 1, 0),
            w(w_)
  {
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
              freq(),
              s_lcp_T(1,0),
              ilist(pars.p.size()),
              ilist_s(pars.p.size()+ 1, 0),
              w(w_)
  {
    // Uploading the frequency file
    std::string tmp_filename = filename + std::string(".occ");
    read_file(tmp_filename.c_str(), freq);
    freq.insert(freq.begin(), 1);

    // Compute the length of the string;
    compute_n();

    // b_p(pfp.n,0);
    verbose("Computing b_p");
    _elapsed_time(compute_b_p());

    verbose("Computing s_lcp_T and ilist");
    _elapsed_time(compute_s_lcp_T_and_ilist());

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

  // Return the frequency of the phrase
  size_t get_freq(size_t phrase){
    return select_ilist_s(phrase + 2) - select_ilist_s(phrase + 1);
  }

  void compute_s_lcp_T_and_ilist(){
    size_t n_phrases = dict.n_phrases();
    size_t len_P = pars.p.size();


    // Check if it is worth to compute it here or read it from file, pushing the computation one step further in the pipeline.

    {
      size_t j = 0;
      ilist_s[j++] = 1; // this is equivalent to set pf.freq[0]=1;
      ilist_s[j] = 1;
      for (size_t i = 1; i < freq.size(); ++i)
      {
        j += freq[i];
        ilist_s[j] = 1;
        freq[i] = 0;
      }
    }

    select_ilist_s = sdsl::bit_vector::select_1_type(&ilist_s);

    size_t last_begin = 0;
    // For all elements of lcpP, compute the corresponding LCP value in T
    for (size_t i = 1; i < len_P; ++i)
    {
      size_t s_i = select_b_p(pars.saP[i] + 1);
      size_t e_i = select_b_p(pars.saP[i] + pars.lcpP[i] + 1);

      size_t l_com_phrases = e_i - s_i;

      assert(pars.lcpP[i] > 0 || l_com_phrases == 0);

      uint_t a = pars.p[pars.saP[i] + pars.lcpP[i]];
      uint_t b = pars.p[pars.saP[i - 1] + pars.lcpP[i]];

      if (a == 0 || b == 0)
      { // we are comparing a phrase with the termination character of P that is not a phrase.
        s_lcp_T.push_back(0);
      }
      else
      {
        // Compute the lcp between phrases a and b
        auto a_in_sa = dict.isaD[dict.select_b_d(a)]; // position of the phrase a in saD
        auto b_in_sa = dict.isaD[dict.select_b_d(b)]; // position of the phrase b in saD

        auto lcp_left = std::min(a_in_sa, b_in_sa) + 1;
        auto lcp_right = max(a_in_sa, b_in_sa);

        size_t lcp_a_b_i = dict.rmq_lcp_D(lcp_left, lcp_right);
        auto lcp_a_b = dict.lcpD[lcp_a_b_i];

        s_lcp_T.push_back(l_com_phrases + lcp_a_b);
      }


      // Update first and last occurrence of each phrase in BWT_P
      size_t prec_phrase_index = (pars.saP[i] == 0 ? len_P : pars.saP[i]) - 1;
      uint_t prec_phrase = pars.p[prec_phrase_index];

      size_t ilist_p = select_ilist_s(prec_phrase + 1) + freq[prec_phrase]++;
      ilist[ilist_p] = i;
    }

    // Computes the ilist for the first element.
    {
      size_t i = 0;
      size_t prec_phrase_index = (pars.saP[i] == 0 ? len_P : pars.saP[i]) - 1;
      uint_t prec_phrase = pars.p[prec_phrase_index];

      size_t ilist_p = select_ilist_s(prec_phrase + 1) + freq[prec_phrase]++;
      ilist[ilist_p] = i;
    }

    rmq_s_lcp_T = sdsl::rmq_succinct_sct<>(&s_lcp_T);
  }

  void clear_unnecessary_elements(){
    // Reducing memory tentative
    pars.isaP.clear(); 
    pars.isaP.shrink_to_fit();
    
    pars.lcpP.clear();
    pars.lcpP.shrink_to_fit();
    
    pars.p.clear();
    pars.p.shrink_to_fit();
    
    freq.clear();
    freq.shrink_to_fit();
  }

  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += dict.serialize(out, child, "dictionary");
    written_bytes += pars.serialize(out, child, "parse");
    // written_bytes += my_serialize(freq, out, child, "frequencies");
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
    // my_load(freq, in);
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
