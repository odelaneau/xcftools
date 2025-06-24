/*******************************************************************************
 * Copyright (C) 2023-2025 Simone Rubinacci
 * Copyright (C) 2023-2025 Olivier Delaneau
 *
 * MIT Licence
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _GTCHECK_H
#define _GTCHECK_H

#include <utils/otools.h>
#include <utils/xcf.h>
#include <utils/bitvector.h>
#include <objects/sparse_genotype.h>
#include <gtcheck/gtcheck_argument_set.h>

struct gtdata
{
	uint32_t nsamples;
	gtcount pop_count;
    int32_t* full_int_buf;
    bitvector binary_bit_buf;
	std::vector<int32_t> sparse_int_buf;

	std::vector<int8_t> unphased_gt;

	int nref;
	int nalt;
	int nhom0;
	int nhom1;
	int nhet;
	int an;
	bool is_phased;

	gtdata() : nsamples(0), full_int_buf(nullptr),
	           nref(0), nalt(0), nhom0(0), nhom1(0), nhet(0), an(0),is_phased(false)
	{
		pop_count.reset();

	}

	gtdata(uint32_t _nsamples) : nsamples(_nsamples), full_int_buf(nullptr),
	           nref(0), nalt(0), nhom0(0), nhom1(0), nhet(0), an(0),is_phased(false)
	{
		pop_count.reset();
		binary_bit_buf.allocate(2 * nsamples);
		sparse_int_buf.resize(2 * nsamples,0);
		unphased_gt.resize(nsamples, 0); // -1 means missing
		full_int_buf = (int32_t*) malloc(2 * nsamples * sizeof(int32_t));
	}

	void initialize(uint32_t _nsamples)
	{
		nsamples = _nsamples;
		pop_count.reset();
		binary_bit_buf.allocate(2 * nsamples);
		sparse_int_buf.resize(2 * nsamples,0);
		unphased_gt.resize(nsamples, 0); // -1 means missing
		if (!full_int_buf)
			full_int_buf = (int32_t*) malloc(2 * nsamples * sizeof(int32_t));
		nref = 0;
		nalt = 0;
		nhom0 = 0;
		nhom1 = 0;
		nhet = 0;
		an = 0;
		is_phased = false;
	}

	void reset_record()
	{
		pop_count.reset();
		sparse_int_buf.clear();
		nref = 0;
		nalt = 0;
		nhom0 = 0;
		nhom1 = 0;
		nhet = 0;
		an = 0;
		is_phased = false;
	}

	void set_sparse(const bool major)
	{
		assert(nsamples >= pop_count.ns + pop_count.mis);
		pop_count.nhom[major] += (nsamples-pop_count.ns-pop_count.mis);
		pop_count.ns = nsamples - pop_count.mis;
	}

	void set_missing()
	{
		pop_count.mis++;
	}
	void set_counts(const bool a0,const bool a1)
	{
		/*
		if (a0==a1) pop_count.nhom[a0] += 2;
		else
		{
			++pop_count.nhet[a0];
			++pop_count.nhet[a1];
		}
		++pop_count.ns;
		*/
		if (a0 == a1)
		    pop_count.nhom[a0] ++;
		else
		    pop_count.nhet ++;

		++pop_count.ns;
	}

	void set_remaining_counts()
	{
		/*
		nref = pop_count.nhet[0] + pop_count.nhom[0];
		nalt = pop_count.nhet[1] + pop_count.nhom[1];
		nhom0 = pop_count.nhom[0];
		nhom1 = pop_count.nhom[1];
		nhet = pop_count.nhet[1];
		an = nref + nalt;
		*/
		nhom0 = pop_count.nhom[0];
		nhom1 = pop_count.nhom[1];
		nhet  = pop_count.nhet;
		nref  = 2*pop_count.nhom[0] + pop_count.nhet;
		nalt  = 2*pop_count.nhom[1] + pop_count.nhet;
		an    = nref + nalt;
	}

	~gtdata()
	{
		if (full_int_buf) free(full_int_buf);
		full_int_buf = nullptr;
		sparse_int_buf.clear();
		unphased_gt.clear();
	}
};
/*
struct ComparisonStats {
	uint64_t n_total = 0;
	uint64_t n_equal = 0;
	uint64_t n_mismatch = 0;
	std::ofstream mismatch_log;

	void start_log(const std::string& filename) {
		mismatch_log.open(filename);
		if (!mismatch_log.is_open()) throw std::runtime_error("Failed to open mismatch log file.");
		mismatch_log << "CHROM\tPOS\tINFO\n";
	}

	void log_mismatch(const std::string& chrom, int pos, const std::string& info) {
		if (mismatch_log.is_open())
			mismatch_log << chrom << "\t" << pos << "\t" << info << "\n";
	}

	void report(std::ostream& out) {
		out << "== Genotype Comparison Report ==\n";
		out << "Total variants compared : " << n_total << "\n";
		out << "Matching records        : " << n_equal << "\n";
		out << "Mismatching records     : " << n_mismatch << "\n";
		out << "Matching %              : " << (100.0 * n_equal / n_total) << "%\n";
	}
};
*/

struct SiteDiff {
    std::string chr, ref, alt;
    int pos;
    std::vector<std::string> diff_fields; // e.g. {"nalt", "nhet"}
    int confusion[3][3]{};
    bool has_confusion = false;
};

struct ComparisonStats {
    uint64_t n_total = 0, n_equal = 0, n_mismatch = 0;

    void add_diff(const bool has_diff)
    {
    	++n_total;
    	has_diff ? ++n_mismatch : ++n_equal;
    }
    void report() {
        vrb.print("== Genotype Comparison Report ==");
        vrb.print("Total variants compared : " + stb.str(n_total));
        vrb.print("Matching records        : " + stb.str(n_equal));
        vrb.print("Mismatching records     : " + stb.str(n_mismatch));
        if (n_total)
        	vrb.print("Matching %              : " + stb.str((100.0 * n_equal / n_total), 2) + "%");
    }
};

class gtcheck {
public:
	const gtcheck_argument_set A;

	std::array<gtdata,2> f;
	//CONSTRUCTOR
	gtcheck(std::vector < std::string > &);
	~gtcheck();

	//ROUTINES
	bool isBCF(std::string);
	bool isXCF(std::string);

	//METHODS
	void run();
	void run_algorithm();

	//
	void read_files_and_initialise();
	void parse_genotypes(xcf_reader& XR, const uint32_t idx_file);
	bool has_gt_difference(const xcf_reader& XR, xcf_writer& XW, bcf1_t* out_rec);
	bool compare_records_exact(xcf_reader &XR, const uint32_t idx0, const uint32_t idx1);
	void hdr_append(bcf_hdr_t* out_hdr);
	void prepare_output(xcf_reader& XR, xcf_writer& XW, const uint32_t idx_file);


	void view(std::vector < std::string > &);
	void write_files_and_finalise();
};

#endif


