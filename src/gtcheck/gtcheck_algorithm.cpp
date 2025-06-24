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

#include <filesystem>
#include <gtcheck/gtcheck_header.h>

void gtcheck::run_algorithm()
{
	tac.clock();
	vrb.title("[GTcheck] Checking XCF files");

	vrb.print("Inititialization");
	vrb.bullet("Opening input files");
	xcf_reader XR(A.mNumThreads, true);
	for (int32_t i=0; i<2; ++i)
	{
		const uint32_t idx_file = XR.addFile(A.mInputFilenames[i]);
		const int32_t typef = XR.typeFile(idx_file);
		vrb.print("  * Opening file [" + A.mInputFilenames[i] + "] (type = " + stb.str(typef) + ")");
		//if (typef != FILE_BINARY) vrb.error("[" + A.mInputFilenames[i] + "] is not a XCF file");
		if (XR.ind_names[idx_file].size() == 0) vrb.error("[" + A.mInputFilenames[i] + "] has no samples");
		f[i].initialize(XR.ind_names[idx_file].size());
	}

	vrb.print("  * Checking sample names");
	if (XR.ind_names[0] != XR.ind_names[1])
	{
		vrb.error("Sample names in the two XCF files do not match: " + A.mInputFilenames[0] + " vs. " + A.mInputFilenames[1]);
	}

	tac.clock();
	vrb.print("Opening input files [done] (" + stb.str(tac.rel_time()/1000) + "s)");

	vrb.bullet("Opening output file");
	xcf_writer XW(A.mOutputFilename, false, A.mNumThreads, false);
	prepare_output(XR,XW,0);
	tac.clock();
	vrb.bullet("Opening output file [done] (" + stb.str(tac.rel_time()/1000) + "s)");

	bcf1_t* out_rec = bcf_init1();
	ComparisonStats stats;
	uint32_t n_variants_total = 0;

	vrb.bullet("Parsing genotypes and checking differences");
	while (XR.nextRecord())
	{
		if ( XR.hasRecord(0) && XR.hasRecord(1) )
		{
			parse_genotypes(XR, 0);
			parse_genotypes(XR, 1);
			const bool has_diff = has_gt_difference(XR, XW, out_rec);
			stats.add_diff(has_diff);
		}
		if (++n_variants_total % 100000 == 0) vrb.bullet("Number of XCF records processed: N = " + stb.str(n_variants_total));
	}
	bcf_destroy(out_rec);
	tac.clock();
	vrb.bullet("Parsing genotypes and checking differences [done] (" + stb.str(tac.rel_time()/1000) + "s)");

	vrb.bullet("Number of variants processed in both files: N = " + stb.str(n_variants_total) +
	           " (shared: " + stb.str(stats.n_total) + ")");

	stats.report();

	XR.close();
	XW.close();
}

void gtcheck::parse_genotypes(xcf_reader& XR, const uint32_t idx_file)
{
	f[idx_file].reset_record();
	const int32_t type = XR.typeRecord(idx_file);
	if (type == RECORD_BCFVCF_GENOTYPE)
	{
	    const int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char** >(&f[idx_file].full_int_buf));
	    f[idx_file].is_phased= false;
		bool missing_seen = false;
	    for (uint32_t i = 0; i < f[idx_file].nsamples; ++i)
	    {
	        const int g0 = f[idx_file].full_int_buf[2*i+0];
	        const int g1 = f[idx_file].full_int_buf[2*i+1];
	        const bool missing = (g0 == bcf_gt_missing || g1 == bcf_gt_missing);
	        if (missing) {
	        	if (!missing_seen) missing_seen = true;
	        	if (f[idx_file].is_phased)
	                vrb.error("Missing data in phased data is not permitted!");
	        	f[idx_file].set_missing();
	            if (A.mDeepCheck)
	                f[idx_file].unphased_gt[i] = -1;
	        } else {
	            if (!missing_seen)
	                f[idx_file].is_phased = (g0 & 1) || (g1 & 1);
	            const bool a0 = (bcf_gt_allele(g0) == 1);
	            const bool a1 = (bcf_gt_allele(g1) == 1);
	            f[idx_file].set_counts(a0, a1);
	            if (A.mDeepCheck)
	                f[idx_file].unphased_gt[i] = a0 + a1;
	        }
	    }
	}
	else if (type == RECORD_BINARY_GENOTYPE) {
		f[idx_file].is_phased=false;
		const int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char* > (&f[idx_file].binary_bit_buf.bytes[0]));
		for(uint32_t i = 0 ; i < f[idx_file].nsamples ; i++)
		{
			const bool a0 = f[idx_file].binary_bit_buf.get(2*i+0);
			const bool a1 = f[idx_file].binary_bit_buf.get(2*i+1);
			const bool missing = (a0 == true && a1 == false);
			missing? f[idx_file].set_missing() : f[idx_file].set_counts(a0, a1);
			if (A.mDeepCheck) f[idx_file].unphased_gt[i] = missing? -1 : (a0 + a1);
		}
	}
	//Convert from binary haplotypes
	else if (type == RECORD_BINARY_HAPLOTYPE)
	{
		f[idx_file].is_phased=true;
		const int32_t n_elements = XR.readRecord(idx_file, reinterpret_cast< char* > (&f[idx_file].binary_bit_buf.bytes[0]));
		for(uint32_t i = 0 ; i < f[idx_file].nsamples ; i++)
		{
			const bool a0 = f[idx_file].binary_bit_buf.get(2*i+0);
			const bool a1 = f[idx_file].binary_bit_buf.get(2*i+1);
			f[idx_file].set_counts(a0,a1);
			if (A.mDeepCheck) f[idx_file].unphased_gt[i] = (a0 + a1);
			//no missing possible? otherwise is_half
		}
	}
	//Convert from sparse genotypes
	else if (type == RECORD_SPARSE_GENOTYPE) {
		f[idx_file].is_phased=false;
		f[idx_file].sparse_int_buf.resize(XR.bin_size[idx_file]/ sizeof(int32_t));
		XR.readRecord(idx_file, reinterpret_cast< char* > (f[idx_file].sparse_int_buf.data()));
		const bool major = (XR.getAF(idx_file)>0.5f);
		if (A.mDeepCheck) std::fill(f[idx_file].unphased_gt.begin(), f[idx_file].unphased_gt.end(), major ? 2 : 0);
		for(uint32_t r = 0 ; r < f[idx_file].sparse_int_buf.size() ; r++)
		{
			sparse_genotype rg(f[idx_file].sparse_int_buf[r]);
			rg.mis ? f[idx_file].set_missing() : f[idx_file].set_counts(rg.al0,rg.al1);
			if (A.mDeepCheck) f[idx_file].unphased_gt[rg.idx] = rg.mis? -1 : (rg.al0 + rg.al1);
		}
		f[idx_file].set_sparse(major);
	}
	else if (type == RECORD_SPARSE_HAPLOTYPE)
	{
		f[idx_file].is_phased=true;
		f[idx_file].sparse_int_buf.resize(XR.bin_size[idx_file]/ sizeof(int32_t));
		XR.readRecord(idx_file, reinterpret_cast< char* > (f[idx_file].sparse_int_buf.data()));
		const bool major = (XR.getAF()>0.5f);
		if (A.mDeepCheck) std::fill(f[idx_file].unphased_gt.begin(), f[idx_file].unphased_gt.end(), major ? 2 : 0);
		for(uint32_t r = 0 ; r < f[idx_file].sparse_int_buf.size() ; r++)
		{
			const int32_t hap_idx = f[idx_file].sparse_int_buf[r];
			const int32_t ind_idx = hap_idx/2;
			const bool a0 = !major;
			const bool a1=(hap_idx%2==0 && r<f[idx_file].sparse_int_buf.size()-1 && f[idx_file].sparse_int_buf[r+1]==hap_idx+1)? a0 : major;//exploits order in sparse vector
			f[idx_file].set_counts( a0, a1);
			if (A.mDeepCheck) f[idx_file].unphased_gt[ind_idx] = (a0 + a1);
			if (a1==a0) ++r;
		}
		f[idx_file].set_sparse( major);
	}
	//Unknown record type
	else vrb.warning("Unrecognized genotype record type [" + stb.str(type) + "] at " + XR.chr + ":" + stb.str(XR.pos));
	f[idx_file].set_remaining_counts();
}
/*
bool gtcheck::compare_records_exact(xcf_reader &XR, const uint32_t idx0, const uint32_t idx1) {
	const int32_t type = XR.typeRecord(idx0);
	if (type != XR.typeRecord(idx1)) return false;

	if (type == RECORD_BCFVCF_GENOTYPE || type == RECORD_BINARY_GENOTYPE || type == RECORD_BINARY_HAPLOTYPE) {
		const int32_t n0 = XR.readRecord(idx0, reinterpret_cast<char*>(&f[0].binary_bit_buf.bytes[0]));
		const int32_t n1 = XR.readRecord(idx1, reinterpret_cast<char*>(&f[1].binary_bit_buf.bytes[0]));
		if (n0 != n1) return false;
		return std::equal(f[0].binary_bit_buf.bytes, f[0].binary_bit_buf.bytes + n0,
		           f[1].binary_bit_buf.bytes);
	}

	if (type == RECORD_SPARSE_GENOTYPE || type == RECORD_SPARSE_HAPLOTYPE) {
		const size_t n_bytes0 = XR.bin_size[idx0];
		const size_t n_bytes1 = XR.bin_size[idx1];
		if (n_bytes0 != n_bytes1) return false;
		f[0].sparse_int_buf.resize(n_bytes0 / sizeof(int32_t));
		f[1].sparse_int_buf.resize(n_bytes1 / sizeof(int32_t));
		XR.readRecord(idx0, reinterpret_cast<char*>(f[0].sparse_int_buf.data()));
		XR.readRecord(idx1, reinterpret_cast<char*>(f[1].sparse_int_buf.data()));
		return f[0].sparse_int_buf == f[1].sparse_int_buf;
	}

	// Unknown or unsupported type
	vrb.warning("Exact record comparison not implemented for type [" + stb.str(type) + "]");
	return false;
}
*/

bool gtcheck::has_gt_difference(const xcf_reader& XR, xcf_writer& XW, bcf1_t* out_rec)
{
	std::vector<std::string> diff_fields;
	//if (valid_record_combinations.count({XR.typeRecord(0), XR.typeRecord(1)}))
	if (f[0].is_phased==f[1].is_phased) // both phased or both unphased
	{
		if (f[0].an          != f[1].an)
			diff_fields.push_back("AN");
		if (f[0].nalt        != f[1].nalt)
			diff_fields.push_back("AC");
		if (f[0].pop_count.mis != f[1].pop_count.mis)
			diff_fields.push_back("NMISS");
		if (f[0].nhom0       != f[1].nhom0)
			diff_fields.push_back("NHOM0");
		if (f[0].nhet        != f[1].nhet)
			diff_fields.push_back("NHET");
		if (f[0].nhom1       != f[1].nhom1)
			diff_fields.push_back("NHOM1");
	}
	else
	{
		int gt_idx = -1, hap_idx = -1;
		if (f[0].is_phased)
		{
			hap_idx = 0;gt_idx = 1;
		}
		else
		{
			hap_idx = 1;gt_idx = 0;
		}
		int diff_nalt  = f[hap_idx].nalt  - f[gt_idx].nalt;
		if (diff_nalt < 0 || diff_nalt > f[gt_idx].pop_count.mis*2)
		    diff_fields.push_back("AC");

		int diff_nhom0 = f[hap_idx].nhom0 - f[gt_idx].nhom0;
		if (diff_nhom0 < 0 || diff_nhom0 > f[gt_idx].pop_count.mis)
		    diff_fields.push_back("NHOM0");

		int diff_nhet  = f[hap_idx].nhet  - f[gt_idx].nhet;
		if (diff_nhet < 0 || diff_nhet > f[gt_idx].pop_count.mis)
		    diff_fields.push_back("NHET");

		int diff_nhom1 = f[hap_idx].nhom1 - f[gt_idx].nhom1;
		if (diff_nhom1 < 0 || diff_nhom1 > f[gt_idx].pop_count.mis)
		    diff_fields.push_back("NHOM1");
	}

    if (A.mDeepCheck)
	{
    	for (size_t i = 0; i < f[0].unphased_gt.size(); ++i) {
    		if (f[0].unphased_gt[i] < 0 || f[1].unphased_gt[i] < 0)
    			continue; // skip missing
            if (f[0].unphased_gt[i] != f[1].unphased_gt[i])
            {
            	if (f[0].unphased_gt[i] != f[1].unphased_gt[i])
            	{
            		diff_fields.push_back("MISMATCH_GT(" + XR.ind_names[0][i] + ")");
            		break;
            	}
            }
        }
    }

	if (diff_fields.empty())
		return false; // Nothing to report

	// Something is different, let's write the record.
	out_rec->rid = bcf_hdr_name2id(XW.hts_hdr, XR.chr.c_str());
	out_rec->pos = XR.pos;
	bcf_update_alleles_str(XW.hts_hdr, out_rec, std::string(XR.ref + "," + XR.alt).c_str());

	std::string field_diff_str = diff_fields.empty() ? "NA" : std::accumulate(std::next(diff_fields.begin()), diff_fields.end(), diff_fields[0],
		[](const std::string &a, const std::string &b) { return a + "," + b; });

	if (bcf_update_info_string(XW.hts_hdr, out_rec, "FD", field_diff_str.c_str()) < 0) {
		vrb.error("Failed to update FD");
		return false;
	}

	int val;
	val = f[0].an;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "AN_F1", &val, 1) < 0) {
		vrb.error("Failed to update AN_F1");
		return false;
	}

	val = f[1].an;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "AN_F2", &val, 1) < 0) {
		vrb.error("Failed to update AN_F2");
		return false;
	}

	val = f[0].nalt;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "AC_F1", &val, 1) < 0) {
		vrb.error("Failed to update AC_F1");
		return false;
	}

	val = f[1].nalt;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "AC_F2", &val, 1) < 0) {
		vrb.error("Failed to update AC_F2");
		return false;
	}

	val = f[0].pop_count.mis;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NMISS_F1", &val, 1) < 0) {
		vrb.error("Failed to update NMISS_F1");
		return false;
	}
	val = f[1].pop_count.mis;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NMISS_F2", &val, 1) < 0) {
		vrb.error("Failed to update NMISS_F2");
		return false;
	}

	val = f[0].nhom0;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NHOMREF_F1", &val, 1) < 0) {
		vrb.error("Failed to update NHOMREF_F1");
		return false;
	}

	val = f[1].nhom0;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NHOMREF_F2", &val, 1) < 0) {
		vrb.error("Failed to update NHOMREF_F2");
		return false;
	}
	val = f[0].nhet;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NHET_F1", &val, 1) < 0) {
		vrb.error("Failed to update NHET_F1");
		return false;
	}
	val = f[1].nhet;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NHET_F2", &val, 1) < 0) {
		vrb.error("Failed to update NHET_F2");
		return false;
	}
	val = f[0].nhom1;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NHOMALT_F1", &val, 1) < 0) {
		vrb.error("Failed to update NHOMALT_F1");
		return false;
	}
	val = f[1].nhom1;
	if (bcf_update_info_int32(XW.hts_hdr, out_rec, "NHOMALT_F2", &val, 1) < 0) {
		vrb.error("Failed to update NHOMALT_F2");
		return false;
	}

	XW.writeRecord(out_rec);
	return true;
}
