/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
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

#ifndef _FILL_TAGS_H
#define _FILL_TAGS_H

#include <utils/otools.h>
#include "fill_tags_argument_set.h"
#include <utils/xcf.h>
#include <containers/bitvector.h>
#include <objects/sparse_genotype.h>

struct AlleleCount {
    std::array<int,2> nhet;
    std::array<int,2> nhom;
    int ns=0;
    int mis=0;

    void reset()
    {
    	nhet={0,0};
    	nhom={0,0};
    	ns=0;
    	mis=0;
    }
};

struct MendelTrio
{
	std::array<int,3> id;
	std::array<int8_t,3> gt;

	MendelTrio(int _id)
	{
		id={_id,-1,-1};
		gt={0,-1,-1};
	}

	void set_gt(const int _id, const int8_t _gt)
	{
		if (_gt > 2 || _gt < -1)
			vrb.error("GT cannot be >2 or <-1");

		for (auto i=0; i<3;++i)
		{
			if (id[i]==_id)
			{
				gt[i]=_gt;
				break;
			}
		}
	}
	void reset(const int8_t _gt)
	{
		if (_gt != 0 && _gt != 2)
			vrb.error("GT cannot be !=0 or 2 in reset");
		for (auto i=0; i<3;++i)
		{
			if (id[i]>=0) gt[i] = _gt;
		}
	}

	int checkMendelError()
	{
		const int8_t kg=gt[0];
		const int8_t fg=gt[1];
		const int8_t mg=gt[2];
		int error = 0;
		if (kg>=0 && fg>=0 && mg>=0) {
			if (fg == 0 && mg == 0 && kg == 1) { error = 1;}
			if (fg == 0 && mg == 0 && kg == 2) { error = 1;}
			if (fg == 0 && mg == 1 && kg == 2) { error = 1;}
			if (fg == 0 && mg == 2 && kg == 0) { error = 1;}
			if (fg == 0 && mg == 2 && kg == 2) { error = 1;}
			if (fg == 1 && mg == 0 && kg == 2) { error = 1;}
			if (fg == 1 && mg == 2 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 0 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 0 && kg == 2) { error = 1;}
			if (fg == 2 && mg == 1 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 2 && kg == 0) { error = 1;}
			if (fg == 2 && mg == 2 && kg == 1) { error = 1;}
		}
		if (kg>=0 && fg>=0 && mg<0) {
			if (fg == 0 && kg == 2) { error = 1;}
			if (fg == 2 && kg == 0) { error = 1;}
		}
		if (kg>=0 && fg<0 && mg>=0) {
			if (mg == 0 && kg == 2) { error = 1;}
			if (mg == 2 && kg == 0) { error = 1;}
		}
		return error;
	}
	int checkMendelTotal(const bool major)
	{
		const int8_t kg=gt[0];
		const int8_t fg=gt[1];
		const int8_t mg=gt[2];
		int total = 0;
		if (kg>=0 && fg>=0 && mg>=0) {
			//total = (!major && kg!=0 && fg!=0 && mg!=0) || (major && kg!=2 && fg!=2 && mg!=2);
			if (major) total = (kg!=2) || (fg!=2) || (mg!=2);
			else total = (kg!=0) || (fg!=0) || (mg!=0);
		}
		if (kg>=0 && fg>=0 && mg<0) {
			//total = (!major && kg!=0 && fg!=0) || (major && kg!=2 && fg!=2);
			if (major) total = (kg!=2) || (fg!=2);
			else total = (kg!=0) || (fg!=0);
		}
		if (kg>=0 && fg<0 && mg>=0) {
			//total = (!major && kg!=0 && mg!=0) || (major && kg!=2 && mg!=2);
			if (major) total = (kg!=2) || (mg!=2);
			else total = (kg!=0) || (mg!=0);
		}
		return total;
	}
};

class fill_tags {
public:
	const fill_tags_argument_set A;

	uint32_t nsamples;

	//oops
    std::vector<std::string> pop_names;
    std::vector<AlleleCount> pop_counts;
    std::vector<std::vector<uint32_t>> pop2samples;//all samples of a population
    std::vector<std::vector<int>> samples2pop;//all populations for each samples

    //mendel
    std::vector<MendelTrio> fam_trio;
    std::vector<std::vector<size_t>> samples2fam;

    //std::vector < std::string > kids;
	//std::vector < int > mothers_idx;
	//std::vector < int > fathers_idx;
	std::vector < int > mendel_errors;
	std::vector < int > mendel_totals;

	//Buffer for input/output

	bitvector binary_bit_buf;
	std::vector<int32_t> sparse_int_buf;

	//CONSTRUCTOR
	fill_tags(std::vector < std::string > &);
	~fill_tags();

	//ROUTINES
	bool isBCF(std::string);
	bool isXCF(std::string);

	//METHODS
	void run();
	void run_algorithm();

	//
	void set_sparse(const uint32_t pop, const bool major);
	void set_missing(const uint32_t pop);
	void set_counts(const uint32_t pop,const bool a0, const bool a1);
	void read_files_and_initialise();
	void hdr_append(bcf_hdr_t* out_hdr);
	void prepare_output(const xcf_reader& XR, xcf_writer& XW,const uint32_t idx_file);
	void process_populations(const xcf_reader& XR, const uint32_t idx_file);
	void parse_genotypes(xcf_reader& XR, const uint32_t idx_file);
	void process_tags(const xcf_reader& XR, xcf_writer& XW,const uint32_t idx_file, std::vector<double>& hwe_probs);
	void calc_hwe(int nref, int nalt, int nhet, std::vector<double>& hwe_probs, float *p_hwe, float *p_exc_het) const;
	void calc_hwe_chisq(const int an, const int fcnt0, const int nhom0, const int nhom1, const int nhet, float *p_chi_square_pval) const;
	void calc_inbreeding_f(const int an, const int fcnt0, const int nhet, float *inbreeding_f) const;
	void calc_mendel_err(int& n_err, int& n_tot, const bool major);

	void process_families(xcf_reader& XR, const uint32_t idx_file);
	void finalize_tags(xcf_reader& XR, const uint32_t idx_file);

	void view(std::vector < std::string > &);
	void write_files_and_finalise();
};

#endif


