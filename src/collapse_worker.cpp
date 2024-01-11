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
#include "collapse_worker.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <sstream>

#include "align_aux.hpp"
#include "aux.hpp"
#include "sequence.hpp"

using namespace gwsc;

void CollapseWorker::operator()() {
    BamBuffer::rpair range;
    while(buffer_->get_next(range)){
        process_range(range);
    }
}

void CollapseWorker::process_range(BamBuffer::rpair range, unsigned int fno){
    fno_filter_ = fno;
    for(auto it = range.first; it != range.second; it++){
        BamDetail & d = *(*it);
        if(fno_filter_ != std::numeric_limits<unsigned int>::max() && fno_filter_ != d.filenum ){
            continue;
        }
        d.xt = bam_aux2A(bam_aux_get(d.b, "RE"));
        if(bhash_){
            tmp_ = bam_aux2Z(bam_aux_get(d.b, "CB"));
            d.barcode = bhash_(tmp_);
        }
        tmp_ = bam_aux2Z(bam_aux_get(d.b, "UB"));
        seq2int<gwsc::ADNA4, uint32_t>(tmp_, d.umi);
        d.make_hash();
        d.lft = d.b->core.pos;
        d.rgt = bam_endpos(d.b) - 1;
        d.processed = false;

        if(bam_is_rev(d.b)){
            d.pos = d.rgt;
        }else{
            d.pos = d.lft;
        }
    }
    sort_reads_index(range.first, range.second, sidx_);
    size_t istart = 0;
    uint32_t lgid = (*(sidx_.front() + range.first))->gid;
    for(size_t i = 0; i < sidx_.size(); i++){
        auto it = range.first + sidx_[i];
        BamDetail & d = *(*it);
        if(d.gid != lgid){
            process_gene_(istart, i, range.first);
            istart = i;
            lgid = d.gid;
        }
    }
    if(istart < sidx_.size()){
        process_gene_(istart, sidx_.size(), range.first);
    }
}

void CollapseWorker::process_range_ds(BamBuffer::rpair range, double ds){
    fds_filter = ds;
    for(auto it = range.first; it != range.second; it++){
        BamDetail & d = *(*it);
        if(fds_filter <= 1.0 && d.rnd >= fds_filter){
            continue;
        }
        d.xt = bam_aux2A(bam_aux_get(d.b, "RE"));
        if(bhash_){
            tmp_ = bam_aux2Z(bam_aux_get(d.b, "CB"));
            d.barcode = bhash_(tmp_);
        }
        tmp_ = bam_aux2Z(bam_aux_get(d.b, "UB"));
        seq2int<gwsc::ADNA4, uint32_t>(tmp_, d.umi);
        d.make_hash();
        d.lft = d.b->core.pos;
        d.rgt = bam_endpos(d.b) - 1;
        d.processed = false;

        if(bam_is_rev(d.b)){
            d.pos = d.rgt;
        }else{
            d.pos = d.lft;
        }
    }
    sort_reads_index(range.first, range.second, sidx_);
    size_t istart = 0;
    uint32_t lgid = (*(sidx_.front() + range.first))->gid;
    for(size_t i = 0; i < sidx_.size(); i++){
        auto it = range.first + sidx_[i];
        BamDetail & d = *(*it);
        if(d.gid != lgid){
            process_gene_(istart, i, range.first);
            istart = i;
            lgid = d.gid;
        }
    }
    if(istart < sidx_.size()){
        process_gene_(istart, sidx_.size(), range.first);
    }
}

void CollapseWorker::process_gene_(size_t start, size_t end, BamBuffer::rit it){
    //std::cout << "    Processing gid = " << gid << " range = " << start << " - " << end << " ffno = " << fno_filter_ << "\n";

    // Correct the UMIs
    size_t rcnt = 0;
    for(size_t i = start; i < end; i++){
        BamDetail & d = *(*(it + sidx_[i]));
        if(fds_filter <= 1.0 && d.rnd >= fds_filter){
            continue;
        }else if(fno_filter_ != std::numeric_limits<unsigned int>::max() && fno_filter_ != d.filenum ){
            continue;
        }
        rcnt++;
    }
    //std::cout << "    rcount = " << rcnt << "\n";
    if(rcnt == 0) return;

    //std::cout << " umap size = " << uhash_.size() << " tot = " << tot << " cor = " << cor << "\n";
    AlignSummary::bint lbarcode = 0;
    bool lbarcode_set = false;
    barcodes_.clear();
    size_t i = start;
    while(i < end){
        BamDetail & d = *(*(it + sidx_[i]));
        if(fds_filter <= 1.0 && d.rnd >= fds_filter){
            i++;
            continue;
        }else if(fno_filter_ != std::numeric_limits<unsigned int>::max() && fno_filter_ != d.filenum ){
            i++;
            continue;
        }
        if(!lbarcode_set) {
            lbarcode_set = true;
            lbarcode = d.barcode;
        }
        if(d.barcode != lbarcode){
            process_barcode_();
            barcodes_.clear();
            lbarcode = d.barcode;
        }
        barcodes_.push_back(&d);
        i++;
    }
    process_barcode_();
}

void CollapseWorker::process_barcode_(){
    if(barcodes_.empty()) return;
    /*
    if(resort){
        sort(barcodes_.begin(), barcodes_.end(),
            [](const BamDetail * p1, const BamDetail * p2) {
                return std::tie(p1->hash, p1->pos) < std::tie(p2->hash, p2->pos);
            });
    }
    */

    //std::cout << "  Barcode " << barcodes_.front()->barcode << "\n";
    umis_.clear();
    uint32_t lumi = barcodes_.front()->umi;
    for(auto d : barcodes_){
        //std::cout << "    " << d->barcode << " " << d->umi << " " << d->pos << "\n";
        if(d->umi != lumi){
            collapse_umi_();
            umis_.clear();
            lumi = d->umi;
        }
        umis_.push_back(d);
    }
    if(!umis_.empty()) collapse_umi_();
}


void CollapseWorker::collapse_umi_(){
    total++;
    rreads += umis_.size();
    build_contigs_();
    fdups_ = 0;
    for(auto & u : umis_){
        if((u->b->core.flag & BAM_FDUP) > 0){
            fdups_++;
        }
    }

    if(build_islands_()){
        //std::cout << "Island group\n";
        fcigar_.clear();
        fquals_.clear();
        fbases_.clear();
        freads_ = 0;
        fsplices_.clear();
        fcoverage_.clear();
        for(size_t i = 0; i < icount_; i++){
            //std::cout << "  Processing group " << i << "\n";
            if(i > 0){
                int d = islands_[i]->lft - islands_[i - 1]->rgt - 1;
                if(d > 0) fcigar_.push_back(CigarElement(d, Cigar::REF_SKIP));
            }
            islands_[i]->merge(umis_, fbases_, fquals_, fcigar_, fsplices_, fcoverage_);
            //islands_[i]->debug(); // __asm__("int $3");
        }

        std::sort(fsplices_.begin(), fsplices_.end());
        fsplices_.erase(std::unique(fsplices_.begin(), fsplices_.end()), fsplices_.end());

        freads_ = umis_.size();
        rdups += fdups_;
        rcollapsed += freads_;
        //std::cout << "freads_ = " << freads_ << " umis size = " << umis_.size() << "\n";
        creads++;

        if((ccount + 1) >= collapsed.size()){
            cbuffer_.get(collapsed, 100);
        }

        bam1_t * bf = umis_[islands_[0]->contigs[0]->index]->b;
        auto & ref = genome_[bf->core.tid];
        unsigned int lft = islands_[0]->lft;
        unsigned int qlft = 0;
        //std::string rstr, qstr, astr;
        unsigned int NM = 0;
        unsigned int bases = 0;
        auto rit = fsplices_.end();
        cgaps_.clear();
        CollapseSplice sp;
        for(auto c : fcigar_){
            switch(c.op){
                case Cigar::MATCH:
                    bases += c.len;
                    for(size_t i = 0; i < c.len; i++, lft++, qlft++){
                        NM += (ref.seq[lft] != fbases_[qlft]);
                        //rstr.push_back(ref.seq[lft]);
                        //qstr.push_back(fbases_[qlft]);
                        //astr.push_back(fbases_[qlft] == ref.seq[lft] ? '|' : '*');
                    }
                    break;
                case Cigar::DEL:
                    /*
                    for(size_t i = 0; i < c.len; i++, lft++){
                        rstr.push_back(ref.seq[lft]);
                        qstr.push_back('_');
                        astr.push_back(' ');
                    }
                    */
                    NM += c.len;
                    lft += c.len;
                    break;
                case Cigar::INS:
                    /*
                    for(size_t i = 0; i < c.len; i++, qlft++){
                        rstr.push_back('_');
                        qstr.push_back(fbases_[qlft]);
                        astr.push_back(' ');
                    }
                    */
                    bases += c.len;
                    NM += c.len;
                    qlft += c.len;
                    break;
                case Cigar::REF_SKIP:
                    sp = {lft - 1, lft + c.len};
                    rit = std::find(fsplices_.begin(), fsplices_.end(), sp);
                    cgaps_.push_back(rit == fsplices_.end() ? '1' : '0');
                    /*
                    rstr += "   ";
                    qstr += "<->";
                    astr += "   ";
                    */
                    lft += c.len;
                    break;
                default:
                    break;
            }
        }

        //ldata.push_back({freads_, bases});

        /*
        //std::cout << "  Fin:  " << fbases_ << "\n";
        std::cout << "  Query: " << qstr << "\n";
        std::cout << "  Align: " << astr << "\n";
        std::cout << "  Ref:   " << rstr << "\n";
        std::cout << "  Qual: ";
        for(size_t i = 0; i < fquals_.size(); i++){
            std::cout << static_cast<char>(fquals_[i] + 33);
        }
        std::cout <<  "\n";
        std::cout << "  Cig : " << fcigar_ << " Reads = " << freads_ << " dups = " << fdups_ << "\n";
        std::cout << "  Qual size = " << fquals_.size() << " seq size = " << fbases_.size() << " cov size = " << fcoverage_.size() << "\n";
        std::cout << "  Cov : ";
        for(auto c : fcoverage_) std::cout << c << ", ";
        std::cout << "\n\n";
        */

        BamDetail * out = collapsed[ccount++];
        make_read_(*out, NM);
    }else{
        rdups += fdups_;
        ambig++;
        rlost += umis_.size();
    }
}

/*
void CollapseWorker::infer_qpos_cigar(BamDetail & d, std::stringstream & ss) {
    d.processed = true;
    int rpos = islands_[0]->lft, qpos = 0;
    bool done = false;
    int clft = 0, crgt = 0;
    for(auto it = fcigar_.begin(); (it != fcigar_.end() && !done); it++){
        const auto & c = *it;
        switch(c.op){
            case Cigar::MATCH:
                for(size_t i = 0; i < c.len; i++){
                    if(rpos == d.lft){
                        clft = qpos;
                    }else if(rpos == d.rgt){
                        crgt = qpos;
                        done = true;
                        break;
                    }
                    rpos++;
                    qpos++;
                }
                break;
            case Cigar::DEL:
                for(size_t i = 0; i < c.len; i++){
                    if(rpos == d.lft){
                        clft = qpos;
                    }else if(rpos == d.rgt){
                        crgt = qpos;
                        done = true;
                        break;
                    }
                    rpos++;
                }
                break;
            case Cigar::INS:
                qpos += c.len;
                break;
            case Cigar::REF_SKIP:
                rpos += c.len;
                break;
            default:
                break;
        }
    }

    if(done){
        ss << clft << "," << crgt;
    }else{
        std::cout << "Couldn't find qpos for collapsed read\n";
    }
}
*/

void CollapseWorker::compress_RLE(const std::vector<uint32_t>& data, std::vector<uint32_t>& compressed){
  //done right before this is called
  // compressed.clear();
  if (data.empty()) {
    return;
  }

  uint32_t prev = data[0];
  uint32_t count = 1;

  for (size_t i = 1; i < data.size(); ++i) {
    if (data[i] == prev) {
      ++count;
    } else {
      compressed.push_back(prev);
      compressed.push_back(count);
      prev = data[i];
      count = 1;
    }
  }

  // Append the last pair
  compressed.push_back(prev);
  compressed.push_back(count);
}

void CollapseWorker::make_read_(BamDetail & bamd, uint32_t NM) {
    bam1_t * bam = bamd.b;
    BamDetail & d = *umis_[islands_[0]->contigs[0]->index];
    
    bamd.filenum = d.filenum;
    bamd.barcode = d.barcode;
    bamd.gid = d.gid;
    bam1_t * bf = d.b;

    bam->core.flag = 0;
    if(bam_is_rev(bf)){
        bam->core.flag |= BAM_FREVERSE;
    }

    bam->core.qual = 255; //a.mapq;  All "uniquely" mapped reads
    bam->core.n_cigar = fcigar_.size();
    bam->core.tid = bf->core.tid;
    bam->core.mtid = -1;
    bam->core.mpos = -1;
    bam->core.isize = -1;
    bam->core.pos = islands_[0]->lft;

    std::stringstream ss;

    fqname_.clear();
    fqname_ += bam_aux2Z(bam_aux_get(bf, "CB"));

    fqname_ += '_';
    fqname_ += std::to_string(umis_[islands_[0]->contigs[0]->index]->gid);
    fqname_ += '_';
    fqname_ += bam_aux2Z(bam_aux_get(bf, "UB"));

    auto qlen = fqname_.size() + 1;
    bam->core.l_extranul = (4 - (qlen & 3)) & 3;
    bam->core.l_qname = qlen + bam->core.l_extranul; // Because of the null terminator
    bam->core.l_qseq  = fbases_.size();

    uint32_t dl = bam->core.n_cigar*4 + bam->core.l_qname + bam->core.l_qseq + (bam->core.l_qseq + 1)/2;

    if(bam->m_data < dl) {
        bam->m_data = dl;
        kroundup32(bam->m_data);
        bam->data = (uint8_t*)realloc(bam->data, bam->m_data);
    }

    bam->l_data = dl;
    // Copy the qname
    uint8_t * p = bam->data;
    memcpy(p, fqname_.c_str(), qlen);
    for(size_t i = 0; i < bam->core.l_extranul; i++){
        p[qlen + i] = 0;
    }
    p += bam->core.l_qname;

    for(auto c : fcigar_){
        uint32_t v = c.packed();
        memcpy(p, &v, 4);
        p+=4;
    }

    for(size_t i = 0; i < fbases_.size(); i++){
        uint8_t * pp = p + (i / 2);
        uint8_t mask = 0xFU << (((~i) & 1) << 2);
        *pp = ((*pp) & ~mask) | seq_nt16_table[static_cast<size_t>(fbases_[i])] << (((~i) & 1) << 2);
    }

    p += (bam->core.l_qseq + 1)/2;
    for(size_t i = 0; i < fquals_.size(); i++){
        *p++ = fquals_[i];
    }

    bam->core.bin = hts_reg2bin(bam->core.pos, bam_endpos(bam), 14, 5);
    if((p - bam->data) != dl){
        std::cout << "Data length wrong dl = " << dl << " compared to " << (p - bam->data) << "\n";
        exit(1);
    }

    /*
    auto ptr = bam_aux_get(bf, "XG");
    if(ptr == NULL) {
        std::cout << "Missing XG tag for read " << bf->data << " freads_ = " << freads_ << "\n";
        //print_tags(bf);
    }
    */
    char * x = bam_aux2Z(bam_aux_get(bf, "XG"));
    bam_aux_append(bam, "XG", 'Z', strlen(x) + 1, reinterpret_cast<uint8_t*>(x));

    x = bam_aux2Z(bam_aux_get(bf, "CB"));
    bam_aux_append(bam, "CB", 'Z', strlen(x) + 1, reinterpret_cast<uint8_t*>(x));

    x = bam_aux2Z(bam_aux_get(bf, "UB"));
    bam_aux_append(bam, "UB", 'Z', strlen(x) + 1, reinterpret_cast<uint8_t*>(x));

    bam_aux_append(bam, "NM", 'i', 4, reinterpret_cast<uint8_t*>(&NM));
    bam_aux_append(bam, "ND", 'i', 4, reinterpret_cast<uint8_t*>(&fdups_));
    bam_aux_append(bam, "NR", 'i', 4, reinterpret_cast<uint8_t*>(&freads_));

    char tmpc = bam_aux2A(bam_aux_get(bf, "RE"));
    bam_aux_append(bam, "RE", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));

    tmpc = bam_aux2A(bam_aux_get(bf, "XS"));
    bam_aux_append(bam, "XS", 'A', 1, reinterpret_cast<uint8_t*>(const_cast<char*>(&tmpc)));

    tmpc = bam_aux2A(bam_aux_get(bf, "RA"));
    bam_aux_append(bam, "RA", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));

    if(!cgaps_.empty()){
        bam_aux_append(bam, "XC", 'Z', cgaps_.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char *>(cgaps_.c_str())));
        //std::cout << "XC = " << cgaps_ << "\n";
    }

    /*
    ss.clear();
    bool hread = false;
    for(size_t i = 0; i < icount_; i++){
        for(auto contig : islands_[i]->contigs){
            if(!umis_[contig->index]->processed){
                if(hread) ss << ';';
                infer_qpos_cigar(*umis_[contig->index], ss);
                hread = true;
            }
        }
    }
    */

    //std::cout << "XR = " << ss.str() << "\n";

    ss.clear();
    bool hread = false;
    for(size_t i = 0; i < icount_; i++){
        auto & isl = *islands_[i];
        for(auto & c : isl.contigs){
            BamDetail & d = *umis_[c->index];
            if(!d.processed){
                bam1_t * bi = d.b;
                if(hread) ss << ';';
                ss << bi->core.pos << ',' << (bam_endpos(bi) - 1) << ',' << c->cig;
                hread = true;
                d.processed = true;
            }
        }
    }

    fqname_ = ss.str();
    bam_aux_append(bam, "XR", 'Z', fqname_.size() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(fqname_.c_str())));




    fcoverageCompressed_.clear();
    compress_RLE(fcoverage_, fcoverageCompressed_);

    bcov_.clear();
    // The base coverage
    // 1 for the array type, 4 for the size, and 4 * fcoverage size for 32 bit ints
    bcov_.resize(5 + fcoverageCompressed_.size() * 4);
    // https://samtools.github.io/hts-specs/SAMv1.pdf
    //I->uint32 t
    bcov_[0] = 'I';
    uint8_t * bp = &bcov_[0];
    (*bp++) = 'I';
    uint32_t l = fcoverageCompressed_.size();
    memcpy(bp, &l, sizeof(uint32_t));
    //bcov_/bp is a char, so need 4 bytes for coverage length
    bp+=4;
    memcpy(bp, &fcoverageCompressed_[0], fcoverageCompressed_.size() * 4);
    bam_aux_append(bam, "CC", 'B', bcov_.size(), reinterpret_cast<uint8_t*>(&bcov_[0]));


    /*
    fqname_.clear();
    for(size_t i = 0; i < fsplices_.size(); i++){
        if(i > 0) fqname_ += '|';
        fqname_ += std::to_string(fsplices_[i].first);
        fqname_ += ',';
        fqname_ += std::to_string(fsplices_[i].second);
    }

    // Splice junctions
    bam_aux_append(bam, "XJ", 'Z', fqname_.size() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(fqname_.c_str())));
    */
}

void CollapseWorker::build_contigs_() {
    ccount_ = 0;
    for(size_t i = 0; i < umis_.size(); i++){
        //unsigned int fs = ccount_;
        ccount_++;
        if(ccount_ >= contigs_.size()) {
            for(size_t i = 0; i < 10; i++) contigs_.push_back(new ReadContig());
        }

        ReadContig * f = contigs_[ccount_-1];
        f->cig.clear();
        f->ins.clear();
        f->lft = umis_[i]->b->core.pos;
        f->rgt = umis_[i]->b->core.pos;
        f->qlft = 0;
        f->index = i;
        f->next = nullptr;
        f->prev = nullptr;

        uint32_t * cig = bam_get_cigar(umis_[i]->b);
        unsigned int cs = 0;
        if(bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP){
            f->qlft += bam_cigar_oplen(cig[0]);
            cs++;
        }
        f->qrgt = f->qlft;
        auto s = bam_get_seq(umis_[i]->b);
        auto q = bam_get_qual(umis_[i]->b);
        for(size_t j = cs; j < umis_[i]->b->core.n_cigar; j++){
            int op =  bam_cigar_op(cig[j]);
            int len = bam_cigar_oplen(cig[j]);
            if(op == BAM_CREF_SKIP){
                // Want one based
                ccount_++;
                if(ccount_ >= contigs_.size()) {
                    for(size_t i = 0; i < 10; i++) contigs_.push_back(new ReadContig());
                }

                auto fn = contigs_[ccount_ - 1];
                fn->cig.clear();
                fn->ins.clear();
                fn->index = i;
                fn->lft = f->rgt + len;
                fn->rgt = fn->lft;
                fn->qlft = f->qrgt;
                fn->qrgt = f->qrgt;
                fn->next = nullptr;
                fn->prev = f;
                f->next = fn;
                f = fn;
            }else if(op == BAM_CMATCH){
                f->qrgt += len;
                f->rgt  += len;
                f->cig.push_back(CigarElement(cig[j]));
            }else if(op == BAM_CDEL){
                f->rgt  += len;
                f->cig.push_back(CigarElement(cig[j]));
            }else if(op == BAM_CINS){
                f->ins.push_back({f->rgt});
                for(int x = 0; x < len; x++, f->qrgt++){
                    f->ins.back().bases.push_back("=ACMGRSVTWYHKDBN"[bam_seqi(s, f->qrgt)]);
                    f->ins.back().quals.push_back(q[f->qrgt]);
                }
                f->cig.push_back(CigarElement(cig[j]));
            }
        }
    }

    // Need to subtract one from each rgt
    for(size_t i = 0; i < ccount_; i++){
        auto & f = *contigs_[i];
        f.rgt--;
        f.qrgt--;
    }

    sort(contigs_.begin(), contigs_.begin() + ccount_,
        [](const ReadContig * p1, const ReadContig * p2) {
            return p1->lft < p2->lft;
    });
}

bool CollapseWorker::build_islands_(){
    icount_ = 1;
    if(icount_ >= islands_.size()) {
        for(size_t i = 0; i < 10; i++) islands_.push_back(new ReadIsland());
    }

    ReadIsland * isl = islands_[icount_ - 1];
    isl->contigs.clear();
    //isl->splices.clear();
    isl->starts.clear();
    isl->ends.clear();
    isl->conns.clear();
    isl->lft = contigs_.front()->lft;
    isl->rgt = contigs_.front()->rgt;
    for(size_t i = 0; i < ccount_; i++){
        auto & f = *contigs_[i];
        //std::cout << "  Check " << isl->lft << " - " << isl->rgt << " vs " << f.lft << " - " << f.rgt
        //    << " c1 = " << (isl->lft <= f.rgt) << " c2 = " << (f.lft <= isl->rgt)
        //    << "\n";
        if(isl->lft <= f.rgt && f.lft <= isl->rgt){
            isl->rgt = std::max(isl->rgt, f.rgt);
            isl->contigs.push_back(contigs_[i]);
            contigs_[i]->island = icount_ - 1;
            if(f.prev != nullptr){
                islands_[f.prev->island]->conns.push_back(f.island);
                //islands_[f.prev->island]->splices.push_back({f.prev->rgt, f.lft});
                islands_[f.prev->island]->ends.push_back(f.prev->rgt);
                islands_[f.island]->starts.push_back(f.lft);
            }
        }else{
            icount_++;
            if(icount_ >= islands_.size()) {
                for(size_t j = 0; j < 10; j++) islands_.push_back(new ReadIsland());
            }
            isl = islands_[icount_ - 1];
            isl->contigs.clear();
            isl->conns.clear();
            //isl->splices.clear();
            isl->starts.clear();
            isl->ends.clear();
            islands_[icount_ - 2]->conns.push_back(icount_ - 1);
            isl->contigs.push_back(contigs_[i]);
            contigs_[i]->island = icount_ - 1;
            isl->lft = f.lft;
            isl->rgt = f.rgt;
            if(f.prev != nullptr){
                //isl->CollapseSplices.push_back({f.prev->lft, f.rgt});
                //islands_[f.prev->island]->splices.push_back({f.prev->rgt, f.lft});
                islands_[f.prev->island]->conns.push_back(f.island);
                islands_[f.prev->island]->ends.push_back(f.prev->rgt);
                islands_[f.island]->starts.push_back(f.lft);
            }
        }
    }

    for(auto i : islands_){
        std::sort(i->conns.begin(), i->conns.end());
        i->conns.erase(std::unique(i->conns.begin(), i->conns.end()), i->conns.end());
        std::sort(i->starts.begin(), i->starts.end());
        i->starts.erase(std::unique(i->starts.begin(), i->starts.end()), i->starts.end());
        std::sort(i->ends.begin(), i->ends.end());
        i->ends.erase(std::unique(i->ends.begin(), i->ends.end()), i->ends.end());
    }
    /*
    for(size_t i = 0; i < umis_.size(); i++){
        std::cout << "  Read: "  << i << " " << debug_read(*umis_[i]) << "\n";
    }
    */
    bool has_skip = false;
    for(size_t i = 0; i < icount_; i++){
        auto & isr = *islands_[i];
        //std::cout << "  Island " << i << " lft " << isr.lft << " - " << isr.rgt << " conns: [";
        bool skip = false;
        if(!isr.conns.empty()){
            //std::cout << isr.conns.front();
            skip |= isr.conns.front() != (i + 1);
            for(size_t j = 1; j < isr.conns.size(); j++){
                //std::cout << ", " << isr.conns[j];
                skip |= isr.conns[i] != (i + 1);
            }
        }

        unsigned int max_sdist = 0, max_edist = 0;

        for(size_t j = 0; j < isr.starts.size(); j++) max_sdist = std::max(isr.starts[j] - isr.lft, max_sdist);
        for(size_t j = 0; j < isr.ends.size(); j++) max_sdist = std::max(isr.rgt - isr.ends[j], max_sdist);
        if(max_sdist > 5 || max_edist > 5){
            skip = true;
        }
        has_skip |= skip;
    }

    return !has_skip;
}
