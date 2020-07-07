/* pfp_thresholds - computing thresholds from prefix free parsing 
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
   \file pfp_thresholds.hpp
   \brief pfp_thresholds.hpp define and build the thresholds from the prefix-free parsing.
   \author Massimiliano Rossi
   \date 02/07/2020
*/

#ifndef _THRESHOLDS_PFP_HH
#define _THRESHOLDS_PFP_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C"
{
#include <gsacak.h>
}

#include <pfp.hpp>

class pfp_thresholds{
public:

    pf_parsing& pf;
    size_t min_s; // Value of the minimum lcp_T in the current run of BWT_T
    size_t pos_s; // Position of the minimum lcp_T in the current run of BWT_T

    uint8_t head; // Head of the current run of BWT_T

    pfp_thresholds(pf_parsing &pfp_, std::string filename) : 
                pf(pfp_),
                min_s(pf.n),
                pos_s(0),
                head(0),
                thresholds(256, pf.n),
                thresholds_pos(256, 0),
                never_seen(256, true)
    {
        // Opening output files
        std::string outfile = filename + std::string(".thr");
        if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".thr_pos");
        if ((thr_pos_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);

        phrase_suffix_t curr;
        phrase_suffix_t prev;

        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {

            if(is_valid(curr)){
                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);  // Store the list of all phrase ids with the same suffix.

                bool same_chars = true;

                phrase_suffix_t next = curr;

                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {

                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_chars = (same_chars && same_suffix.back().bwt_char == next.bwt_char);

                        same_suffix.push_back(next);
                    }

                }

                // Simple case
                if (same_chars)
                {

                    for (auto curr : same_suffix)
                    {
                        // curr = elem;
                        // Compute phrase boundary lcp
                        int_t lcp_suffix = compute_lcp_suffix(curr,prev);

                        // Update min_s
                        update_min_s(lcp_suffix,j);

                        update_bwt(curr.bwt_char, pf.get_freq(curr.phrase));
                        update_min_s(lcp_suffix,j);

                        // get min_s_P using rmq_s_lcp_T
                        size_t left = pf.select_ilist_s(curr.phrase + 1);
                        size_t right = pf.select_ilist_s(curr.phrase + 2) - 1;
                        if (left < right && lcp_suffix >= curr.suffix_length)
                        {
                            size_t pos_s_P = pf.rmq_s_lcp_T(left + 1, right);
                            size_t min_s_P = pf.s_lcp_T[pos_s_P] + curr.suffix_length - pf.dict.length_of_phrase(curr.phrase);
                            // Update current min
                            update_min_s(min_s_P, j + pos_s_P - left);
                        }

                        // Update prevs
                        prev = curr;

                        j += pf.get_freq(curr.phrase);
                    }
                }
                else
                {
                    // Hard case
                    int_t lcp_suffix = compute_lcp_suffix(curr,prev);

                    typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;

                    // using lambda to compare elements.
                    auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                        return *lhs.first > *rhs.first;
                    };

                    std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                    for (auto s: same_suffix)
                    {
                        size_t begin = pf.select_ilist_s(s.phrase + 1);
                        size_t end = pf.select_ilist_s(s.phrase + 2);
                        pq.push({&pf.ilist[begin], {&pf.ilist[end], s.bwt_char}});
                    }

                    size_t prev_occ;
                    bool first = true;
                    while (!pq.empty())
                    {
                        auto curr_occ = pq.top();
                        pq.pop();

                        if (!first)
                        {
                            // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
                            lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
                        }
                        first = false;
                        // Update min_s
                        update_min_s(lcp_suffix, j);
                        
                        update_bwt(curr_occ.second.second, 1);
                        update_min_s(lcp_suffix, j);

                        // Update prevs
                        prev_occ = *curr_occ.first;

                        // Update pq
                        curr_occ.first++;
                        if (curr_occ.first != curr_occ.second.first)
                            pq.push(curr_occ);

                        j += 1;
                    }

                    prev = same_suffix.back();
                }
                curr = next;
                
            }
            else
            {
                inc(curr);
            }
        }

        // Close output files
        fclose(thr_file);
        fclose(thr_pos_file);
    }

private:
    typedef struct
    {
        size_t i = 0; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t phrase = 0;
        size_t suffix_length = 0;
        int_da sn = 0;
        uint8_t bwt_char = 0;
    } phrase_suffix_t;

    
    size_t j = 0;

    std::vector<uint64_t> thresholds;
    std::vector<uint64_t> thresholds_pos;
    std::vector<bool> never_seen;

    FILE *thr_file;
    FILE *thr_pos_file;

    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= pf.dict.saD.size())
            return false;
        s.sn = pf.dict.saD[s.i];
        s.phrase = pf.dict.rank_b_d(s.sn); 
        // s.phrase = pf.dict.daD[s.i] + 1; // + 1 because daD is 0-based
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pf.ilist.size()));
        s.suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(s.sn + 1) + 1) - s.sn - 1;
        if(is_valid(s))
            s.bwt_char = (s.sn == pf.w ? 0 : pf.dict.d[s.sn - 1]);
        return true;
    }

    inline bool is_valid(phrase_suffix_t& s)
    {
        // avoid the extra w # at the beginning of the text
        if (s.sn < pf.w)
            return false;
        // Check if the suffix has length at least w and is not the complete phrase.
        if (pf.dict.b_d[s.sn] != 0 || s.suffix_length < pf.w)
            return false;
        
        return true;
    }
    
    inline int_t min_s_lcp_T(size_t left, size_t right)
    {
        // assume left < right
        if (left > right)
            std::swap(left, right);

        assert(pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] >= pf.w);

        return (pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] - pf.w);
    }

    inline int_t compute_lcp_suffix(phrase_suffix_t& curr, phrase_suffix_t& prev)
    {
        int_t lcp_suffix = 0;

        if (j > 0)
        {
            // Compute phrase boundary lcp
            lcp_suffix = pf.dict.lcpD[curr.i];
            for (size_t k = prev.i + 1; k < curr.i; ++k)
            {
                lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
            }

            if (lcp_suffix >= curr.suffix_length && curr.suffix_length == prev.suffix_length)
            {
                // Compute the minimum s_lcpP of the phrases following the two phrases
                // we take the first occurrence of the phrase in BWT_P
                size_t left = pf.ilist[pf.select_ilist_s(curr.phrase + 1)]; //size_t left = first_P_BWT_P[phrase];
                // and the last occurrence of the previous phrase in BWT_P
                size_t right = pf.ilist[pf.select_ilist_s(prev.phrase + 2) - 1]; //last_P_BWT_P[prev_phrase];
                
                lcp_suffix += min_s_lcp_T(left,right);
            }
        }

        return lcp_suffix;
    }

    inline void update_min_s(int_t val, size_t pos)
    {
        if (val < min_s)
        {
            min_s = val;
            pos_s = j;
        }
    }

    // We can put here the check if we want to store the LCP or stream it out
    inline void new_min_s(int_t val, size_t pos)
    {
        min_s = val;
        pos_s = j;
    }

    inline void update_bwt(uint8_t next_char, size_t length)
    {
        if (head != next_char)
        {
            // Print threshold
            print_threshold(next_char);
            head = next_char;
            // lengths.push_back(0); // Debug only
            // Create the new min
            new_min_s(pf.n + 10, j);
        }

        // lengths.back() += length; // Debug only
    }



    inline void print_threshold(uint8_t next_char){

        // Update thresholds
        for (auto character: pf.dict.alphabet)
        {
            if (character == head) continue;
            if (min_s < thresholds[character])
            {
                thresholds[character] = min_s;
                thresholds_pos[character] = pos_s;
            }
        }

        if (never_seen[next_char])
        {
            // if(next_char >= 1) // I am assuming that the BWT has only one 0
                never_seen[next_char] = false;
        }
        else
        {
            if (fwrite(&thresholds[next_char], THRBYTES, 1, thr_file) != 1)
                error("SA write error 1");
            if (fwrite(&thresholds_pos[next_char], THRBYTES, 1, thr_pos_file) != 1)
                error("SA write error 1");
        }


        thresholds[next_char] = pf.n;
    }

};

#endif /* end of include guard: _THRESHOLDS_PFP_HH */
