/* pfp_thresholds - Test of construction of MEM thresholds from prefix free parsing data structures
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
   \file pfp_thresholds.cpp
   \brief pfp_thresholds.cpp test construction of MEM-thresholds from prefix-free parsing.
   \author Massimiliano Rossi
   \date 17/06/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <pfp.hpp>

#include <malloc_count.h>


int main(int argc, char* const argv[]) {


  Args args;
  parseArgs(argc, argv, args);

  // TODO: Include cmd line option to load from file/compute
  // verbose("Loading PFP data structures from file");
  // std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  // pf_parsing<> pf;
  // std::string filename = args.filename + pf.filesuffix();
  // sdsl::load_from_file(pf, filename);

  // Computing prefix-free parsing
  verbose("Window size set to: " , args.w);

  verbose("Computing PFP data structures");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

  pf_parsing pf(args.filename, args.w);

  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("PFP DS construction complete");
  verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
  // auto time = std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count();


  verbose("Building the thresholds");
  

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

  // This code gets timed

  // Building the sampled LCP array of T in corrispondence of the beginning of each phrase.
  verbose("Building the thresholds - sampled LCP");

  size_t n_phrases = pf.dict.n_phrases();
  size_t len_P = pf.pars.p.size();

  std::vector<int_t>  s_lcp_T(1,0);                  // LCP array of T sampled in corrispondence of the beginning of each phrase.
  // std::vector<uint_t> first_P_BWT_P(n_phrases+1); // Position of the first occurrence of the phrase i in BWT_P
  // std::vector<uint_t> last_P_BWT_P(n_phrases+1);  // Position of the last occurrence of the phrase i in BWT_P
  std::vector<uint_t> min_s_P; // Value of the minimum s_lcp_T in each run of P[SA_P]
  std::vector<uint_t> pos_s_P; // Position of the minimum s_lcp_T in each run of P[SA_P]

  std::vector<std::vector<uint_t>> occs_P_BWT_P(n_phrases+1);           // Positions of the occurrences of the phrase i in BWT_P

  // std::vector<bool> first_P_BWT_P_visit(n_phrases+1,false); // Position of the first occurrence of the phrase i in BWT_P

  size_t last_begin = 0;
  // For all elements of lcpP, compute the corresponding LCP value in T
  for(size_t i = 1; i < len_P; ++i){
    size_t s_i = pf.select_b_p(pf.pars.saP[i] + 1);
    size_t e_i = pf.select_b_p(pf.pars.saP[i] + pf.pars.lcpP[i] + 1);

    size_t l_com_phrases = e_i - s_i;

    assert(pf.pars.lcpP[i]>0 || l_com_phrases==0);

    uint_t a = pf.pars.p[pf.pars.saP[i] + pf.pars.lcpP[i]];
    uint_t b = pf.pars.p[pf.pars.saP[i-1] + pf.pars.lcpP[i]];

    if(a == 0 || b == 0){ // we are comparing a phrase with the termination character of P that is not a phrase.
      s_lcp_T.push_back(0);
    }else{
      // Compute the lcp between phrases a and b
      auto a_in_sa = pf.dict.isaD[pf.dict.select_b_d(a)]; // position of the phrase a in saD
      auto b_in_sa = pf.dict.isaD[pf.dict.select_b_d(b)]; // position of the phrase b in saD

      auto lcp_left = std::min(a_in_sa,b_in_sa) + 1;
      auto lcp_right = max(a_in_sa,b_in_sa);

      size_t lcp_a_b_i = pf.dict.rmq_lcp_D(lcp_left, lcp_right);
      auto lcp_a_b = pf.dict.lcpD[lcp_a_b_i];

      s_lcp_T.push_back(l_com_phrases + lcp_a_b);
    }


    // Compute mins and poss
    if (pf.pars.p[pf.pars.saP[i-1]] == pf.pars.p[pf.pars.saP[i]]){
      // If both suffixes starts with the same phrase, update the value of the min
      if(min_s_P.back() > s_lcp_T.back() - l_com_phrases){
        min_s_P.back() = s_lcp_T.back() - l_com_phrases;
        pos_s_P.back() = i - last_begin;
      }
    } else {
      min_s_P.push_back(pf.n + 10);
      pos_s_P.push_back(0);
      last_begin = i;
    }

    // Update first and last occurrence of each phrase in BWT_P
    size_t prec_phrase_index = (pf.pars.saP[i] == 0? len_P: pf.pars.saP[i]) - 1;
    uint_t prec_phrase = pf.pars.p[prec_phrase_index];

    // last_P_BWT_P[prec_phrase] = i;
    // if(!first_P_BWT_P_visit[prec_phrase]){
    //   first_P_BWT_P_visit[prec_phrase] = true;
    //   first_P_BWT_P[prec_phrase] = i;
    // }
    occs_P_BWT_P[prec_phrase].push_back(i);
  }

  // Computing the thresholds
  verbose("Building the thresholds - min_s and pos_s");
  sdsl::rmq_succinct_sct<> rmq_s_lcp_T = sdsl::rmq_succinct_sct<>(&s_lcp_T);

  std::vector<size_t> min_s(1,pf.n); // Value of the minimum lcp_T in each run of BWT_T
  std::vector<size_t> pos_s(1,0); // Position of the minimum lcp_T in each run of BWT_T
  // std::vector<uint_t> min_s(1,pf.n); // Value of the minimum lcp_T in each run of BWT_T
  // std::vector<uint_t> pos_s(1,0); // Position of the minimum lcp_T in each run of BWT_T

  std::vector<uint8_t> heads(1,0);
  std::vector<size_t> lengths(1,0); //  Debug only

  size_t prev_i = 0;
  size_t prev_phrase = 0;
  size_t prev_suffix_length = 0;

  assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);
  size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
  size_t j = 0;
  while (i < pf.dict.saD.size()){

    auto sn = pf.dict.saD[i];
    // Check if the suffix has length at least w and is not the complete phrase.
    auto phrase = pf.dict.daD[i] + 1; // + 1 because daD is 0-based
    assert(phrase > 0 && phrase < pf.freq.size()); 
    size_t suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(sn + 1) + 1) - sn - 1;

    if (sn < pf.w){ // avoid the extra w # at the beginning of the text
      i = i + 1;
    }else if (pf.dict.b_d[sn] == 0 && suffix_length >= pf.w){

      // Compute the next character of the BWT of T
      uint8_t bwt_char = (sn == pf.w ? 0 : pf.dict.d[sn - 1]);
      std::vector<size_t> same_suffix(1, phrase);           // Store the list of all phrase ids with the same suffix.
      std::vector<uint8_t> bwt_chars(1, bwt_char); // the character corresponding to the phrase id
      bool same_chars = true;

      size_t next = i + 1;
      while (next < pf.dict.saD.size() && (pf.dict.lcpD[next] >= suffix_length)){
        auto next_sn = pf.dict.saD[next];
        // Check if the suffix has length at least w and is not the complete phrase.
        auto next_phrase = pf.dict.daD[next] + 1; // + 1 because daD is 0-based
        assert(next_phrase > 0 && next_phrase < pf.freq.size());
        size_t next_suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(next_sn + 1) + 1) - next_sn - 1;

        assert(next_suffix_length >= suffix_length);
        assert((pf.dict.b_d[next_sn] == 0 && next_suffix_length >= pf.w) || (next_suffix_length != suffix_length));
        if (next_suffix_length == suffix_length){
          same_suffix.push_back(next_phrase);
          bwt_char = (next_sn == pf.w ? 0 : pf.dict.d[next_sn - 1]);
          same_chars = (same_chars && bwt_chars.back() == bwt_char);
          bwt_chars.push_back( bwt_char );
        }

        next = next + 1;
      }


      // Simple case
      if(same_chars){
        
        for( auto phrase: same_suffix){

          int_t lcp_suffix = 0;

          if (j > 0)
          {
            // Compute phrase boundary lcp
            lcp_suffix = pf.dict.lcpD[i];
            for (size_t k = prev_i + 1; k < i; ++k)
            {
              lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
            }

            if (lcp_suffix >= suffix_length && suffix_length == prev_suffix_length)
            {
              // Compute the minimum s_lcpP of the phrases following the two phrases
              // we take the first occurrence of the phrase in BWT_P
              size_t left = occs_P_BWT_P[phrase][0]; //size_t left = first_P_BWT_P[phrase];
              // and the last occurrence of the previous phrase in BWT_P
              size_t right = occs_P_BWT_P[prev_phrase].back(); //last_P_BWT_P[prev_phrase];
              // assume left < right
              if (left > right)
                std::swap(left, right);

              assert(s_lcp_T[rmq_s_lcp_T(left + 1, right)] >= pf.w);

              lcp_suffix += s_lcp_T[rmq_s_lcp_T(left + 1, right)] - pf.w;
            }
          }

          // Update min_s
          if (lcp_suffix < min_s.back())
          {
            min_s.back() = lcp_suffix;
            pos_s.back() = j;
          }

          auto next_BWT_char = bwt_chars[0];

          if (heads.back() != next_BWT_char)
          {
            heads.push_back(next_BWT_char);
            lengths.push_back(0); // Debug only
            // Create the new min
            min_s.push_back(lcp_suffix);
            pos_s.push_back(j);
          }

          lengths.back() += pf.freq[phrase]; // Debug only

          // Update current min
          if (lcp_suffix >= suffix_length && min_s.back() > suffix_length + min_s_P[phrase])
          {
            min_s.back() = min_s_P[phrase + suffix_length];
            pos_s.back() = pos_s_P[phrase];
          }

          // Update prevs
          prev_i = i;
          prev_phrase = phrase;
          prev_suffix_length = suffix_length;

          j += pf.freq[phrase];
          i += 1;
        }
      }else{
        // Hard case
        int_t lcp_suffix = 0;

        if (j > 0)
        {
          // Compute phrase boundary lcp
          lcp_suffix = pf.dict.lcpD[i];
          for (size_t k = prev_i + 1; k < i; ++k)
          {
            lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
          }

          if (lcp_suffix >= suffix_length && suffix_length == prev_suffix_length)
          {
            // Compute the minimum s_lcpP of the phrases following the two phrases
            // we take the first occurrence of the phrase in BWT_P
            size_t left = occs_P_BWT_P[phrase][0]; //size_t left = first_P_BWT_P[phrase];
            // and the last occurrence of the previous phrase in BWT_P
            size_t right = occs_P_BWT_P[prev_phrase].back(); //last_P_BWT_P[prev_phrase];
            // assume left < right
            if (left > right)
              std::swap(left, right);

            assert(s_lcp_T[rmq_s_lcp_T(left + 1, right)] >= pf.w);

            lcp_suffix += s_lcp_T[rmq_s_lcp_T(left + 1, right)] - pf.w;
          }
        }

        typedef std::pair < std::vector<uint_t>::iterator, std::pair<std::vector<uint_t>::iterator, uint8_t> > pq_t;

        // using lambda to compare elements.
        auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
          return *lhs.first > *rhs.first;
        };

        std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
        for(size_t k = 0; k < same_suffix.size(); ++k){
          auto phrase = same_suffix[k];
          pq.push({occs_P_BWT_P[phrase].begin(), {occs_P_BWT_P[phrase].end(),bwt_chars[k]}});
        }

        size_t prev_occ;
        bool first = true;
        while(!pq.empty()){
          auto curr = pq.top(); pq.pop();

          if(!first){
            // Compute the minimum s_lcpP of the phrases following the two phrases
            // we take the current occurrence of the phrase in BWT_P
            size_t left = *curr.first;
            // and the previous occurrence of the previous phrase in BWT_P
            size_t right = prev_occ;
            // assume left < right
            if (left > right)
              std::swap(left, right);

            assert(s_lcp_T[rmq_s_lcp_T(left + 1, right)] >= pf.w);

            lcp_suffix = suffix_length + s_lcp_T[rmq_s_lcp_T(left + 1, right)] - pf.w;       
          }
          first = false;
          // Update min_s
          if (lcp_suffix < min_s.back())
          {
            min_s.back() = lcp_suffix;
            pos_s.back() = j;
          }

          auto next_BWT_char = curr.second.second;

          if (heads.back() != next_BWT_char)
          {
            heads.push_back(next_BWT_char);
            lengths.push_back(0); // Debug only
            // Create the new min
            min_s.push_back(lcp_suffix);
            pos_s.push_back(j);
          }

          lengths.back()++; // Debug only

          // Update prevs
          prev_occ = *curr.first;

          // Update pq
          curr.first ++;
          if(curr.first != curr.second.first)
            pq.push(curr);

          j += 1;
        }

        i = next;
        prev_i = i - 1;
        prev_phrase = same_suffix.back();
        prev_suffix_length = suffix_length;
      }
      i = next;
      assert(i==next);
    } else {
      i = i+1;
    }

  }

  verbose("Building the thresholds - constructing thresholds");

  // Opening output files
  FILE *thr_file;
  std::string outfile = args.filename + std::string(".thr");
  if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");
  
  FILE *thr_pos_file;
  outfile = args.filename + std::string(".thr_pos");
  if ((thr_pos_file = fopen(outfile.c_str(), "w")) == nullptr)
    error("open() file " + outfile + " failed");

  // std::vector<size_t> thresholds;
  // std::vector<size_t> thresholds_pos_s;
  std::vector<uint64_t> last_seen(256, 0);
  std::vector<bool> never_seen(256, true);

  sdsl::rmq_succinct_sct<> rmq_min_s = sdsl::rmq_succinct_sct<>(&min_s);


  for(size_t i = 1; i < heads.size(); ++i){
    if(never_seen[heads[i]]){
      never_seen[heads[i]] = false;
    }else{
      size_t j = rmq_min_s(last_seen[heads[i]]+1, i-1);
      if (fwrite(&min_s[j], THRBYTES, 1, thr_file) != 1)
        error("SA write error 1");
      if (fwrite(&pos_s[j], THRBYTES, 1, thr_pos_file) != 1)
        error("SA write error 1");
      // thresholds.push_back(min_s[j]);
      // thresholds_pos_s.push_back(pos_s[j]);
    }
    last_seen[heads[i]] = i;
  }

  // Close output files
  fclose(thr_file);
  fclose(thr_pos_file);

  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration<double, std::ratio<1>>(t_end - t_start).count();
  verbose("Elapsed time (s): ", time);

  auto mem_peak = malloc_count_peak();
  verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo)
  {
    // space = thresholds.size() * sizeof(thresholds[0]);
    // space += thresholds_pos_s.size() * sizeof(thresholds_pos_s[0]);
    verbose("Thresholds size (bytes): ", space);
  }

  // verbose("Storing the Thresholds to file");
  // std::string outfile = args.filename + std::string(".thr");
  // write_file(outfile.c_str(), thresholds);
  // outfile = args.filename + std::string(".thr_pos");
  // write_file(outfile.c_str(), thresholds_pos_s);
  
  if (args.store)
  {
    verbose("Storing the BWT to file");
    std::vector<uint8_t> bwt;
    for(int i = 1; i < heads.size(); ++i){
      bwt.insert(bwt.end(), lengths[i], heads[i] );
    }
    outfile = args.filename + std::string(".my.bwt");
    write_file(outfile.c_str(), bwt);
  }

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
  }
