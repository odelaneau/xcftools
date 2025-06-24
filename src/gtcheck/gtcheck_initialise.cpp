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
#include <utils/otools.h>
#include <utils/xcf.h>
#include <utils/bitvector.h>
#include <objects/sparse_genotype.h>

void gtcheck::read_files_and_initialise()
{
	rng.setSeed(A.mSeed);
}

void gtcheck::hdr_append(bcf_hdr_t* out_hdr)
{
	bcf_hdr_printf(out_hdr, "##INFO=<ID=FD,Number=1,Type=String,Description=\"Differing fields in the two files\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=AN_F1,Number=1,Type=Integer,Description=\"AN in file 1\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=AN_F2,Number=1,Type=Integer,Description=\"AN in file 2\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=AC_F1,Number=1,Type=Integer,Description=\"AC in file 1\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=AC_F2,Number=1,Type=Integer,Description=\"AC in file 2\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NMISS_F1,Number=1,Type=Integer,Description=\"NMISS in file 1\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NMISS_F2,Number=1,Type=Integer,Description=\"NMISS in file 2\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NHOMREF_F1,Number=1,Type=Integer,Description=\"NHOM REF in file 1\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NHOMREF_F2,Number=1,Type=Integer,Description=\"NHOM REF in file 2\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NHET_F1,Number=1,Type=Integer,Description=\"NHET in file 1\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NHET_F2,Number=1,Type=Integer,Description=\"NHET in file 2\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NHOMALT_F1,Number=1,Type=Integer,Description=\"NHOM ALT in file 1\">");
	bcf_hdr_printf(out_hdr, "##INFO=<ID=NHOMALT_F2,Number=1,Type=Integer,Description=\"NHOM ALT in file 2\">");
	//bcf_hdr_printf(out_hdr, "##FORMAT=<ID=GT_F1,Number=1,Type=Integer,Description=\"Genotype field in file 1\">");
	//bcf_hdr_printf(out_hdr, "##FORMAT=<ID=GT_F2,Number=1,Type=Integer,Description=\"Genotype field in file 2\">");
}

//FIXME XR should be const here.
void gtcheck::prepare_output(xcf_reader& XR, xcf_writer& XW, const uint32_t idx_file)
{
    // 1. Start a clean header and set fileformat first
    bcf_hdr_t* out_hdr = bcf_hdr_init("w");
    bcf_hdr_append(out_hdr, "##fileformat=VCFv4.3");
    bcf_hdr_append(out_hdr, "##source=XCFtools_gtcheck");
    // 2. Copy only contig lines from input
    bcf_hdr_t* in_hdr = XR.sync_reader->readers[idx_file].header;
    for (int i = 0; i < in_hdr->n[BCF_DT_CTG]; ++i) {
        const char* contig_name = bcf_hdr_id2name(in_hdr, i);
        int contig_len = in_hdr->id[BCF_DT_CTG][i].val->info[0];
        std::ostringstream oss;
        if (contig_len > 0)
            oss << "##contig=<ID=" << contig_name << ",length=" << contig_len << ">";
        else
            oss << "##contig=<ID=" << contig_name << ">";
        bcf_hdr_append(out_hdr, oss.str().c_str());
    }
    // 3. Add custom INFO and FORMAT fields
    hdr_append(out_hdr);
    // 4. Duplicate for XW and flush()
    XW.hts_hdr = bcf_hdr_dup(out_hdr);
    bcf_hdr_add_sample(XW.hts_hdr, NULL);
    bcf_hdr_sync(XW.hts_hdr);
    // 5. Write header
    if (bcf_hdr_write(XW.hts_fd, XW.hts_hdr) < 0)
        helper_tools::error("Failing to write BCF/header");
    if (!XW.hts_fidx.empty())
        if (bcf_idx_init(XW.hts_fd, XW.hts_hdr, 14, XW.hts_fidx.c_str()))
            helper_tools::error("Initializing .csi");
    bcf_clear1(XW.hts_record);

    bcf_hdr_destroy(out_hdr);
}

