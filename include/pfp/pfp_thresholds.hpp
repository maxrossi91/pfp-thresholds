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
    size_t  length = 0; // Length of the current run of BWT_T

    bool rle;

    pfp_thresholds(pf_parsing &pfp_, std::string filename, bool rle_ = false) : 
                pf(pfp_),
                min_s(pf.n),
                pos_s(0),
                head(0),
                sa_mod(pf.n - pf.w + 1ULL),
                thresholds(256, pf.n),
                thresholds_pos(256, 0),
                never_seen(256, true),
                rle(rle_)
    {
        // Opening output files
        std::string outfile = filename + std::string(".thr");
        if ((thr_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".thr_pos");
        if ((thr_pos_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".ssa");
        if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".esa");
        if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        if(rle)
        {
            outfile = filename + std::string(".bwt.heads");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

            outfile = filename + std::string(".bwt.len");
            if ((bwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

        }else{

            outfile = filename + std::string(".bwt");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

        }

        outfile = filename + std::string(".slcp");
        if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
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

                init_first_last_occ();
                update_first_last_occ(curr);

                phrase_suffix_t next = curr;

                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {

                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_chars = (same_chars && same_suffix.back().bwt_char == next.bwt_char);

                        same_suffix.push_back(next);

                        update_first_last_occ(next);
                    }

                }

                // Simple case
                if (same_chars)
                {

                    update_ssa(same_suffix[0], first_occ);

                    for (auto curr : same_suffix)
                    {
                        // curr = elem;
                        // Compute phrase boundary lcp
                        int_t lcp_suffix = compute_lcp_suffix(curr,prev);

                        // Update min_s
                        update_min_s(lcp_suffix,j);
                        update_lcp(lcp_suffix);

                        update_bwt(curr.bwt_char, pf.get_freq(curr.phrase));
                        update_min_s(lcp_suffix,j);


                        // Update prevs
                        prev = curr;

                        j += pf.get_freq(curr.phrase);
                    }

                    update_esa(same_suffix[0], last_occ);
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
                        size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                        size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                        pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
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
                        update_lcp(lcp_suffix);

                        update_ssa(curr, *curr_occ.first);

                        update_bwt(curr_occ.second.second, 1);
                        update_min_s(lcp_suffix, j);

                        update_esa(curr, *curr_occ.first);
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

        // print lat BWT char and SA sample
        print_sa();
        print_bwt();

        // Close output files
        fclose(thr_file);
        fclose(thr_pos_file);
        fclose(ssa_file);
        fclose(esa_file);
        fclose(bwt_file);
        if(rle)
            fclose(bwt_file_len);
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

    size_t first_occ = pf.n; // First occurrence of same suffix phrases in BWT_P
    size_t last_occ = 0;     // Last occurrence of same suffix phrases in BWT_P

    const size_t sa_mod = (pf.n - pf.w + 1ULL);

    size_t ssa = 0;
    size_t esa = 0;

    size_t prev_lcp = 0;
    size_t curr_lcp = 0;

    std::vector<uint64_t> thresholds;
    std::vector<uint64_t> thresholds_pos;
    std::vector<bool> never_seen;

    FILE *bwt_file;
    FILE *bwt_file_len;

    FILE *ssa_file;
    FILE *esa_file;

    FILE *thr_file;
    FILE *thr_pos_file;

    FILE *lcp_file;

    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= pf.dict.saD.size())
            return false;
        s.sn = pf.dict.saD[s.i];
        s.phrase = pf.dict.rank_b_d(s.sn); 
        // s.phrase = pf.dict.daD[s.i] + 1; // + 1 because daD is 0-based
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pf.pars.ilist.size()));
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
                size_t left = pf.pars.ilist[pf.pars.select_ilist_s(curr.phrase + 1)]; //size_t left = first_P_BWT_P[phrase];
                // and the last occurrence of the previous phrase in BWT_P
                size_t right = pf.pars.ilist[pf.pars.select_ilist_s(prev.phrase + 2) - 1]; //last_P_BWT_P[prev_phrase];
                
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

    inline void update_ssa(phrase_suffix_t &curr, size_t pos)
    {   // We do not need to add w because pf.pos_T has w character more at the beginning
        ssa = (sa_mod + pf.pos_T[pos] - curr.suffix_length) % (sa_mod); // + pf.w;
        // ssa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(ssa < (pf.n - pf.w + 1ULL));
    }

    inline void update_esa(phrase_suffix_t &curr, size_t pos)
    {
        esa = (sa_mod + pf.pos_T[pos] - curr.suffix_length)% (sa_mod);// + pf.w;
        // esa = (pf.pos_T[pos] - curr.suffix_length)% (pf.n - pf.w + 1ULL);// + pf.w;
        assert(esa < (pf.n - pf.w + 1ULL));
    }

    inline void update_lcp(size_t _lcp)
    {
        curr_lcp = _lcp;
    }

    inline void print_sa()
    {
        if (j < (pf.n - pf.w + 1ULL))
        {
            size_t pos = j;
            if (fwrite(&pos, SSABYTES, 1, ssa_file) != 1)
                error("SA write error 1");
            if (fwrite(&ssa, SSABYTES, 1, ssa_file) != 1)
                error("SA write error 2");
        }

        if(j > 0)
        {
            size_t pos = j-1;
            if (fwrite(&pos, SSABYTES, 1, esa_file) != 1)
                error("SA write error 1");
            if (fwrite(&esa, SSABYTES, 1, esa_file) != 1)
                error("SA write error 2");
        }

    }

    inline void print_lcp()
    {
        if (fwrite(&curr_lcp, SSABYTES, 1, lcp_file) != 1)
            error("LCP write error 1");

    }

    inline void print_bwt()
    {
        if(length > 0)
        {
            if(rle)
            {
                // Write the head
                if (fputc(head, bwt_file) == EOF)
                    error("BWT write error 1");
                
                // Write the length
                if (fwrite(&length, BWTBYTES, 1, bwt_file_len) != 1)
                    error("BWT write error 2");
            }else{

                for(size_t i = 0; i < length; ++i)
                {
                    if (fputc(head, bwt_file) == EOF)
                        error("BWT write error 1");
                }
                // uint8_t* run = new uint8_t[length];
                // memset(run,head,length);

                // if (fwrite(run, 1, length, bwt_file) != length)
                //     error("BWT write error 1");

                // delete run;
            }

        }

    }

    inline void update_bwt(uint8_t next_char, size_t length_)
    {
        if (head != next_char)
        {
            // Print threshold
            print_threshold(next_char);
            print_sa();
            print_bwt();
            print_lcp();

            prev_lcp = curr_lcp;
            head = next_char;
            // lengths.push_back(0); // Debug only
            length = 0; // Debug only
            // Create the new min
            new_min_s(pf.n + 10, j);
        }

        length += length_;
        // lengths.back() += length; // Debug only
    }

    inline void init_first_last_occ()
    {
        first_occ = pf.n; // First occurrence of same suffix phrases in BWT_P
        last_occ = 0;     // Last occurrence of same suffix phrases in BWT_P
    }

    inline void update_first_last_occ(phrase_suffix_t &curr) 
    {
        size_t begin = pf.pars.select_ilist_s(curr.phrase + 1);
        size_t end = pf.pars.select_ilist_s(curr.phrase + 2) - 1;

        size_t first = pf.pars.ilist[begin];
        size_t last  = pf.pars.ilist[end];

        if(first_occ > first)
            first_occ = first;

        if(last_occ < last)
            last_occ = last;
    }
    
    
    inline void print_threshold(uint8_t next_char)
    {

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

            // Write a zero so the positions of thresholds and BWT runs are the same
            size_t zero = 0;
            if (fwrite(&zero, THRBYTES, 1, thr_file) != 1)
                error("THR write error 1");
            if (fwrite(&zero, THRBYTES, 1, thr_pos_file) != 1)
                error("THR write error 2");
        }
        else
        {
            if (fwrite(&thresholds[next_char], THRBYTES, 1, thr_file) != 1)
                error("THR write error 3");
            if (fwrite(&thresholds_pos[next_char], THRBYTES, 1, thr_pos_file) != 1)
                error("THR write error 4");
        }


        thresholds[next_char] = pf.n;
    }
};

#endif /* end of include guard: _THRESHOLDS_PFP_HH */
