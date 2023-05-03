#pragma once
/*
Copyright (c) 2018-2020 Gavin W. Wilson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "bam_genes.hpp"
#include "aux.hpp"
#include "align_aux.hpp"

#include <thread>
#include <list>
#include <array>

namespace gwsc {

using PileupSplice = std::pair<unsigned int, unsigned int>;
using PileupSplices = std::vector<PileupSplice>;

struct PileupRead {
    void init(BamDetail & d);

    bool next();

    bool between() const {
        return it == cigar.end();
    }

    uint8_t base() const {
        //     "=ACMGRSVTWYHKDBN";
        return "4014244434444444"[bam_seqi(bam_get_seq(d->b), qpos)] - '0';
    }

    uint8_t qual() const {
        return *(bam_get_qual(d->b) + qpos);
    }

    // Given the 'CC' tag from the modified coverage, which gives the number of
    // raw (uncollapsed) reads associated with the SNV, and the total number.
    // Two optimizations make this more complicated
    // 1) Assumed the total number of reads for each base is less than
    // uint16_t=65535, so packed both (number of positive reads) and (total
    // number of reads) into 1 uint32_t. Here I pass along that uint32_t without
    // splitting
    // 2) Since so many neighboring bases have the same values, I did a
    // rudimentary read length encoding where I store (value, # times repeated).
    // Could do gzip, it's all small number of elements. So here, given a target
    // `qpos`, I iterate through to find the result. Could destructure once and
    // just do a vector look up, dunno. Could also do a gzip decompression

    //TODO don't read in the bam every time, store the umi_coverage_tmp


    uint32_t make_umi_coverage() const {
      unsigned int CC_len = bam_auxB_len(bam_aux_get(d->b, "CC"));
      std::vector<uint32_t> umi_coverage_tmp;
      for(unsigned int CC_i=0; CC_i < CC_len; CC_i++){
        umi_coverage_tmp.push_back(
                                   bam_auxB2i(bam_aux_get(d->b, "CC"),
                                              CC_i));
      }

      //B/c of read length encoding, should always come in pairs 
      if((CC_len % 2) != 0){
            __asm__("int $3");
      }

      auto umi_coverage_it = umi_coverage_tmp.begin();
      unsigned int currentI = 0;
      uint32_t  currentCovg = 0;

      //Really should start currentI at -1, and then add 1 in the first round to
      //get to 0. But b/c of wanting an unsigned int so no casting, start at 0
      //and then need qpos to be +1.
      while(currentI < (qpos+1)){
        currentCovg = *(umi_coverage_it++);
        currentI   += *(umi_coverage_it++);
      }

      return currentCovg;
    }

    int tid() const {
        return d->b->core.tid;
    }

    int qlen() const {
        return d->b->core.l_qseq;
    }


    std::vector<unsigned int>     qpositions;
    CigarString                   cigar;
    CigarString::const_iterator   it;
    BamDetail                   * d;
    unsigned int                  cpos;
    unsigned int                  qpos;
    unsigned int                  rpos;
    unsigned int                  rend;
    unsigned int                  ibases;
    unsigned int                  dbases;
    unsigned int                  abases;
    unsigned int                  NM;
    unsigned int                  qdist;
    char                          strand;
};

struct PileupOut{
    PileupOut(){

    }

    PileupOut(PileupRead & r) 
        : 
            d(r.d), rpos(r.rpos), qpos(r.qpos), qdist(r.qdist),
            base(r.base()), qual(r.qual()), ibases(r.ibases), 
            dbases(r.dbases), abases(r.abases), NM(r.NM), rev(r.strand == '-'),
            umi_coverage(r.make_umi_coverage())
    {

    }

    BamDetail *       d;
    unsigned int      rpos;
    unsigned int      qpos;
    unsigned int      qdist;
    uint8_t           base;
    uint8_t           qual;
    unsigned int      ibases;
    unsigned int      dbases;
    unsigned int      abases;
    unsigned int      NM;
    bool              rev;
    uint32_t          umi_coverage;
};

//template <typename T>
class Pileup {
    public:
        using const_iterator = std::vector<PileupOut>::const_iterator;

        Pileup(bool use_dups = false) 
            : dups_(use_dups) {

        }

        ~Pileup(){
            for(auto & r : reads_){
                delete r;
            }
            reads_.clear();
        }

        void build_reads(BamBuffer::rpair range);
        bool next(std::vector<PileupOut> & out);

        int                tid = -1;
        unsigned int       pos = 0;

    private:
        void get_reads_();
        void pileup_reads_(std::vector<PileupOut> & out);

        std::list<PileupRead*>       buffer_;
        std::vector<PileupRead*>     reads_;
        size_t                       idx_ = 0;
        size_t                       count_ = 0;
        int                          ctid_ = -1;
        unsigned int                 cpos_ = 0;
        bool                         dups_;
};

}
