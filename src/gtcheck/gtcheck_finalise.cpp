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


#include <gtcheck/gtcheck_header.h>

void gtcheck::write_files_and_finalise() {
	vrb.title("Finalization:");

	/*
	const std::string bfname = helper_tools::get_name_from_vcf(A.mOutputFilename) + ".mendel.ind.txt.gz";
	//Mendel per sample summary
	vrb.title("Writing Mendel per sample summary in [" + bfname + "]");
	output_file fds(bfname + ".mendel.ind.txt.gz");
	for (int kidx = 0 ; kidx < XR.ind_names[idx_file].size() ; kidx++)
	{
		fds << 	XR.ind_names[idx_file][kidx] << "\t" <<
				XR.ind_fathers[idx_file][kidx] << "\t" <<
				XR.ind_mothers[idx_file][kidx] << "\t" <<
				XR.ind_pops[idx_file][kidx] << "\t" <<
				mendel_errors[kidx] << "\t" <<
				mendel_totals_fam_all[kidx] << "\t" <<
				mendel_totals_fam_minor[kidx] << std::endl;
	}
	fds.close();
	*/


	//step0: Measure overall running time
	vrb.bullet("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
