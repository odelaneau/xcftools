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

#ifndef _CONCAT_H
#define _CONCAT_H

#include <utils/otools.h>
#include <utils/xcf.h>

class concat {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//FILE DATA
	int nfiles;

	std::vector < std::string > filenames;
	std::vector < int > prev_readers;

	//SAMPLE DATA

	int nsamples;

	std::array<int,2> nswap;
	std::array<std::vector<bool>,2> swap_phase;
	std::vector < int > nmatch;
	std::vector < int > nmism;

	std::vector < int > nsites_buff_d2;

	int32_t *GTa, *GTb, mGTa, mGTb;

	//CONSTRUCTOR
	concat();
	~concat();

	//PARAMETERS
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void concatenate(std::vector < std::string > &);
	void read_files_and_initialise();
	void run();
	void concat_naive();
	void ligate();
	void write_files_and_finalise();
	//Helpers
	void concat_naive_check_headers(xcf_writer& XW, const std::string& fname);
	void check_hrecs(const bcf_hdr_t *hdr0, const bcf_hdr_t *hdr, const char *fname0, const char *fname);
	void scan_overlap(const int ifname,const char* seek_chr, int seek_pos);
	void update_distances();
	void phase_update(bcf_hdr_t *hdr, bcf1_t *line, const bool uphalf);
	void write_record(htsFile *, bcf_hdr_t * ,  bcf_hdr_t * ,bcf1_t *, const bool uphalf);
};

#endif


