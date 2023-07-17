/*******************************************************************************
 * Copyright (C) 2022-2023 Simone Rubinacci
 * Copyright (C) 2022-2023 Olivier Delaneau
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
#include <concat/concat_header.h>
#include <utils/xcf.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>
#include <utils/otools.h>
#include <utils/basic_stats.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define GET(n,i)	(((n)>>(i))&1U)
#define TOG(n,i)	((n)^=(1UL<<(i)))

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

void concat::run() {
	if (options.count("naive"))
	{
		concat_naive();
	}
	else
	{

	}
}

void concat::concat_naive()
{
	tac.clock();
	concat_naive_check_headers();
	vrb.title("Concatenating files:");
	std::string fname = options["output"].as < std::string > ();
	xcf_writer XW(fname, false, 1);
	uint64_t offset_seek = 0;

    for (size_t i=0; i<filenames.size(); i++)
    {
    	tac.clock();
    	vrb.print2("  * Concatenating " + filenames[i]);
    	xcf_reader XR(1);
    	const int32_t idx_file = XR.addFile(filenames[i]);
    	const int32_t type = XR.typeFile(idx_file);
    	if (type != FILE_BINARY) vrb.error("[" + filenames[i] + "] is not a XCF file");

    	if (i==0)
    	{
    		XW.writeHeader(XR.sync_reader->readers[0].header);
    		if (!std::filesystem::exists(stb.remove_extension(filenames[i]) + ".fam")) vrb.error("File does not exists: " + stb.remove_extension(filenames[i]) + ".bin");
    		std::ifstream fam_ifile(stb.remove_extension(filenames[i]) + ".fam");
    		std::ofstream fam_ofile(stb.remove_extension(fname) + ".fam");
    		fam_ofile << fam_ifile.rdbuf();
    		fam_ofile.close();
			fam_ifile.close();
    	}

    	while (XR.nextRecord())
    	{
    		XW.writeInfo(XR.chr, XR.pos, XR.ref, XR.alt, XR.rsid, XR.getAC(), XR.getAN());
    		XW.writeSeekField(XR.bin_type[0], XR.bin_seek[0] + offset_seek, XR.bin_size[0]);
    	}

        if (!std::filesystem::exists(stb.remove_extension(filenames[i]) + ".bin")) vrb.error("File does not exists: " + stb.remove_extension(filenames[i]) + ".bin");
		std::ifstream bin_ifile(stb.remove_extension(filenames[i]) + ".bin", std::ios::in | std::ios::binary);
		XW.bin_fds << bin_ifile.rdbuf();
		XW.bin_fds.flush();
		offset_seek = XW.bin_fds.tellp();
		bin_ifile.close();
        vrb.print("\t(" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
    }
    XW.bin_fds.close();
    XW.close();

    vrb.print("Writing completed.");
}

// This is a C++ friendly modification of vcfconcat.c from bcftools.
// Copyright (C) 2013-2023 Genome Research Ltd.
// Author: Petr Danecek <pd3@sanger.ac.uk>
void concat::concat_naive_check_headers()
{
    vrb.print2("  * Checking the headers of " + stb.str(nfiles) + " files");
    bcf_hdr_t *hdr0 = NULL;
    int i,j;
    for (i=0; i<filenames.size(); i++)
    {
        htsFile *fp = hts_open(filenames[i].c_str(), "r"); if ( !fp ) vrb.error("Failed to open: " + filenames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) vrb.error("Failed to parse header: " + filenames[i]);
        htsFormat type = *hts_get_format(fp);
        hts_close(fp);

        if ( i==0 )
        {
            hdr0 = hdr;
            continue;
        }

        // check the samples
        if ( bcf_hdr_nsamples(hdr0)!=bcf_hdr_nsamples(hdr) )
        	vrb.error("Cannot concatenate, different number of samples: " + stb.str(bcf_hdr_nsamples(hdr0)) + " vs " + stb.str(bcf_hdr_nsamples(hdr0)) + " in "+ filenames[0] + " vs " + filenames[i]);
         for (j=0; j<bcf_hdr_nsamples(hdr0); j++)
            if ( strcmp(hdr0->samples[j],hdr->samples[j]) )
            	vrb.error("Cannot concatenate, different samples in "+ filenames[0] + " vs " + filenames[i]);

        // if BCF, check if tag IDs are consistent in the dictionary of strings
        if ( type.compression!=bgzf )
            vrb.print("The --naive option works only for compressed BCFs as main file for the XCF file format, sorry :-/\n");


        check_hrecs(hdr0,hdr,filenames[0].c_str(),filenames[i].c_str());
        check_hrecs(hdr,hdr0,filenames[i].c_str(),filenames[0].c_str());

        bcf_hdr_destroy(hdr);
    }
    if ( hdr0 ) bcf_hdr_destroy(hdr0);
    vrb.print(". Done, the headers are compatible.");
}

// This is a C++ friendly modification of vcfconcat.c from bcftools.
// Copyright (C) 2013-2023 Genome Research Ltd.
// Author: Petr Danecek <pd3@sanger.ac.uk>
void concat::check_hrecs(const bcf_hdr_t *hdr0, const bcf_hdr_t *hdr, const char *fname0, const char *fname)
{
    int j;
    for (j=0; j<hdr0->nhrec; j++)
    {
        bcf_hrec_t *hrec0 = hdr0->hrec[j];
        if ( hrec0->type!=BCF_HL_FLT && hrec0->type!=BCF_HL_INFO && hrec0->type!=BCF_HL_FMT && hrec0->type!=BCF_HL_CTG ) continue;    // skip fiels w/o IDX
        int itag = bcf_hrec_find_key(hrec0, "ID");
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, hrec0->type, "ID", hrec0->vals[itag], NULL);

        std::string type;
        if ( hrec0->type==BCF_HL_FLT ) type = "FILTER";
        if ( hrec0->type==BCF_HL_INFO ) type = "INFO";
        if ( hrec0->type==BCF_HL_FMT ) type = "FORMAT";
        if ( hrec0->type==BCF_HL_CTG ) type = "contig";

        if ( !hrec )
        	vrb.error("Cannot use --naive, incompatible headers, the tag " + type + "/" + stb.str(hrec0->vals[itag]) + " not present in " + stb.str(fname));

        int idx0 = bcf_hrec_find_key(hrec0, "IDX");
        int idx  = bcf_hrec_find_key(hrec,  "IDX");
        if ( idx0<0 || idx<0 )
        	vrb.error("fixme: unexpected IDX<0 for " + type + "/" + stb.str(hrec0->vals[itag]) + " in " + stb.str(fname0) + " or " + stb.str(fname));
        if ( strcmp(hrec0->vals[idx0],hrec->vals[idx]) )
        	vrb.error("Cannot use --naive. different order the tag in " + type + "/" + stb.str(hrec0->vals[itag]) + " in " + stb.str(fname0) + " vs " + stb.str(fname));
    }
}

// This is a C++ friendly modification of vcfconcat.c from bcftools.
// Copyright (C) 2013-2023 Genome Research Ltd.
// Author: Petr Danecek <pd3@sanger.ac.uk>
void concat::ligate()
{
	tac.clock();
	vrb.title("Ligating chunks");
	std::string fname = options["output"].as < std::string > ();
	xcf_writer XW(fname, false, 1);
	uint64_t offset_seek = 0;

	bcf_srs_t * sr =  bcf_sr_init();
	sr->require_index = 1;
	int n_threads = options["thread"].as < int > ();
	if (n_threads > 1) if (bcf_sr_set_threads(sr, n_threads) < 0) vrb.error("Failed to create threads");

	bcf_hdr_t * out_hdr = NULL;
	bcf1_t *line = bcf_init();
	std::vector<int> start_pos(nfiles);

	for (int f = 0, prev_chrid = -1 ; f < nfiles ; f ++)
	{
		htsFile *fp = hts_open(filenames[f].c_str(), "r"); if ( !fp ) vrb.error("Failed to open: " + filenames[f] + ".");
		bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) vrb.error("Failed to parse header: " + filenames[f] +".");
		out_hdr = bcf_hdr_merge(out_hdr,hdr);
        if ( bcf_hdr_nsamples(hdr) != bcf_hdr_nsamples(out_hdr) ) vrb.error("Different number of samples in " + filenames[f] + ".");
        for (int j=0; j<bcf_hdr_nsamples(hdr); j++)
        	if ( std::string(out_hdr->samples[j]) != std::string(hdr->samples[j]) )  vrb.error("Different sample names in " + filenames[f] + ".");

        int ret = bcf_read(fp, hdr, line);
		if ( ret!=0 ) vrb.error("Empty file detected: " + filenames[f] +".");
        else
        {
            int chrid = bcf_hdr_id2int(out_hdr,BCF_DT_CTG,bcf_seqname(hdr,line));
            start_pos[f] = chrid==prev_chrid ? line->pos : -1;
            prev_chrid = chrid;
        }
        bcf_hdr_destroy(hdr);
        if ( hts_close(fp)!=0 ) vrb.error("Close failed: " + filenames[f] + ".");
	}
    for (int i=1; i<nfiles; i++) if ( start_pos[i-1]!=-1 && start_pos[i]!=-1 && start_pos[i]<start_pos[i-1] ) vrb.error("The files not in ascending order");
    int i = 0, nrm = 0;
	nsamples = bcf_hdr_nsamples(out_hdr);
	nswap = {0,0};
	swap_phase = {std::vector<bool>(nsamples, false), std::vector<bool>(nsamples, false)};
	nmatch = std::vector < int > (nsamples, 0);
	nmism = std::vector < int > (nsamples, 0);

	GTa = GTb = NULL;
	mGTa = 0, mGTb=0;

	XW.writeHeader(out_hdr);
	if (!std::filesystem::exists(stb.remove_extension(filenames[i]) + ".fam")) vrb.error("File does not exists: " + stb.remove_extension(filenames[i]) + "fam");
	std::ifstream fam_ifile(stb.remove_extension(filenames[i]) + ".fam");
	std::ofstream fam_ofile(stb.remove_extension(fname) + ".fam");
	fam_ofile << fam_ifile.rdbuf();
	fam_ofile.close();
	fam_ifile.close();

	int n_variants = 0;
	int n_variants_at_start_cnk = 0;
	line = bcf_init();
	int chunk_counter=0;
	int n_sites_buff = 0;

	int prev_readers_size = 0;
    std::string prev_chr = "";
    std::array<int,2> prev_pos = {0,0};
    int first_pos = 0;
    int ifname = 0;

	vrb.bullet("#samples = " + stb.str(nsamples));
	vrb.print("");
	tac.clock();

    while ( ifname < nfiles )
    {
        int new_file = 0;
        while ( sr->nreaders < 2 && ifname < nfiles )
        {
            if ( !bcf_sr_add_reader (sr, filenames[ifname].c_str())) vrb.error("Failed to open " + filenames[ifname] + ".");
            new_file = 1;
            ifname++;
            if ( start_pos[ifname-1]==-1 ) break;   // new chromosome, start with only one file open
            if ( ifname < nfiles && start_pos[ifname]==-1 ) break; // next file starts on a different chromosome
        }

        // is there a line from the previous run? Seek the newly opened reader to that position
        int seek_pos = -1;
        int seek_chr = -1;
        if ( bcf_sr_has_line(sr,0) )
        {
            bcf1_t *line0 = bcf_sr_get_line(sr,0);
            bcf_sr_seek(sr, bcf_seqname(sr->readers[0].header,line0), line0->pos);
            seek_pos = line0->pos;
            seek_chr = bcf_hdr_name2id(out_hdr, bcf_seqname(sr->readers[0].header,line0));
        }
        else if ( new_file ) bcf_sr_seek(sr,NULL,0);  // set to start

        int nret;
        while ( (nret = bcf_sr_next_line(sr)) )
        {
            if ( !bcf_sr_has_line(sr,0) ) if ( bcf_sr_region_done(sr,0) )  bcf_sr_remove_reader(sr, 0);

            // Get a line to learn about current position
            for (i=0; i<sr->nreaders; i++) if ( bcf_sr_has_line(sr,i) ) break;
            bcf1_t *line = bcf_sr_get_line(sr,i);

            // This can happen after bcf_sr_seek: indel may start before the coordinate which we seek to.
            if ( seek_chr>=0 && seek_pos>line->pos && seek_chr==bcf_hdr_name2id(out_hdr, bcf_seqname(sr->readers[i].header,line)) ) continue;
            seek_pos = seek_chr = -1;

            //  Check if the position overlaps with the next, yet unopened, reader
            int must_seek = 0;
            while ( ifname < nfiles && start_pos[ifname]!=-1 && line->pos >= start_pos[ifname] )
            {
                must_seek = 1;
                if ( !bcf_sr_add_reader(sr, filenames[ifname].c_str())) vrb.error("Failed to open " + filenames[ifname] + ".");
                if  (sr->nreaders>2) vrb.error("Three files overlapping at position: " + std::to_string(line->pos+1));
                ifname++;
            }
            if ( must_seek )
            {
                bcf_sr_seek(sr, bcf_seqname(sr->readers[i].header,line), line->pos);
                seek_pos = line->pos;
                seek_chr = bcf_hdr_name2id(out_hdr, bcf_seqname(sr->readers[i].header,line));
                continue;
            }

            if ( sr->nreaders>1 && ((!bcf_sr_has_line(sr,0) && !bcf_sr_region_done(sr,0)) || (!bcf_sr_has_line(sr,1) && !bcf_sr_region_done(sr,1))) )
            {
            	const bool uphalf = !bcf_sr_has_line(sr,0);
            	line = bcf_sr_get_line(sr,uphalf);
            	//TODO//write_record(out_fp, out_hdr, sr->readers[0].header, line,uphalf);
                prev_pos[uphalf]=line->pos;
                prev_readers_size = sr->nreaders;
				n_variants++;
            	continue;
            }

            if (sr->nreaders<2)
            {
            	if (prev_readers_size == 0)
            	{
            		n_variants_at_start_cnk = n_variants;
            		prev_chr = std::string(bcf_seqname(sr->readers[i].header,line));
            		first_pos = (int)bcf_sr_get_line(sr,i)->pos + 1;
            		vrb.wait("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-]");
            	}
            	else if (prev_readers_size == 2)
    			{
            		n_variants_at_start_cnk = n_variants;
            		prev_chr = std::string(bcf_seqname(sr->readers[i].header,line));
            		first_pos = (int)bcf_sr_get_line(sr,i)->pos + 1;
    				vrb.wait("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-]");
            		//after a buffer we go back to one reader. Chunk 1 is now chunk 0.
    				n_sites_buff = 0;
            		nswap[0]=nswap[1];
            		swap_phase[0] = swap_phase[1];
    			}
				line = bcf_sr_get_line(sr,i);
				//TODO//write_record(out_fp, out_hdr, sr->readers[0].header, line,i);
            	prev_pos[i]=line->pos;
            }
            else
            {
            	if (n_sites_buff==0)
            	{
            		prev_chr = std::string(bcf_seqname(sr->readers[i].header,line));
    				vrb.print("Cnk " + stb.str(ifname-2) + " [" + prev_chr + ":" + stb.str(first_pos) + "-" + stb.str(prev_pos[0] + 1) + "] [L=" + stb.str(n_variants-n_variants_at_start_cnk) + "]" );
            		scan_overlap(ifname, bcf_seqname(sr->readers[i].header,line), line->pos);
            	}

				const bool uphalf = n_sites_buff >= nsites_buff_d2.back();
				line = bcf_sr_get_line(sr,uphalf);
				//TODO//write_record(out_fp, out_hdr, sr->readers[uphalf].header, line, uphalf);
				++n_sites_buff;
	            prev_pos[0]=prev_pos[1]=line->pos;
            }
            prev_readers_size = sr->nreaders;
    		n_variants++;
        }
        if ( sr->nreaders ) while ( sr->nreaders ) bcf_sr_remove_reader(sr, 0);
    }
	vrb.print("Cnk " + stb.str(ifname-1) + " [" + prev_chr + ":" + stb.str(first_pos) + "-" + stb.str(prev_pos[0] + 1) + "] [L=" + stb.str(n_variants-n_variants_at_start_cnk) + "]" );
	bcf_hdr_destroy(out_hdr);
	bcf_sr_destroy(sr);
	if (line) bcf_destroy(line);
	free(GTa);
    free(GTb);

	if (n_variants == 0) vrb.error("No variants to be phased in files");
	XW.bin_fds.close();
	XW.close();

	vrb.title("Writing completed [L=" + stb.str(n_variants) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void concat::phase_update(bcf_hdr_t *hdr, bcf1_t *line, const bool uphalf)
{
	//TODO: next GT HERE: need to do sparse
    int i, nGTs = bcf_get_genotypes(hdr, line, &GTa, &mGTa);
    if ( nGTs <= 0 ) return;    // GT field is not present
    for (i=0; i<bcf_hdr_nsamples(hdr); i++)
    {
		if ( !swap_phase[uphalf][i] ) continue;
		int *gt = &GTa[i*2];
		if ( bcf_gt_is_missing(gt[0]) || gt[1]==bcf_int32_vector_end ) continue;
        if (!bcf_gt_is_phased(gt[0]) || !bcf_gt_is_phased(gt[1])) continue;
        const int gt0 = bcf_gt_phased(bcf_gt_allele(gt[1])==1);
        const int gt1 = bcf_gt_phased(bcf_gt_allele(gt[0])==1);
        gt[0] = gt0;
        gt[1] = gt1;
    }
    bcf_update_genotypes(hdr,line,GTa,nGTs);
}

void concat::update_distances()
{
	//TODO: next GT HERE: needs to be sparse
	for (int i = 0 ; i < nsamples; i++)
	{
	    int *gta = &GTa[i*2];
	    int *gtb = &GTb[i*2];
	    if ( gta[1]==bcf_int32_vector_end || gtb[1]==bcf_int32_vector_end ) continue;
	    if ( bcf_gt_is_missing(gta[0]) || bcf_gt_is_missing(gta[1]) || bcf_gt_is_missing(gtb[0]) || bcf_gt_is_missing(gtb[1]) ) continue;
	    if ( !bcf_gt_is_phased(gta[1]) || !bcf_gt_is_phased(gtb[1]) ) continue;
	    if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gta[1]) || bcf_gt_allele(gtb[0])==bcf_gt_allele(gtb[1]) ) continue;
	    if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[0]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[1]) )
	    {
	        if ( swap_phase[0][i] ) nmism[i]++; else nmatch[i]++;
	    }
	    if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[1]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[0]) )
	    {
	        if ( swap_phase[0][i] ) nmatch[i]++; else nmism[i]++;
	    }
	}
}

void concat::write_record(htsFile *fd, bcf_hdr_t * out_hdr, bcf_hdr_t * hdr_in, bcf1_t *line, const bool uphalf)
{
	//TODO: next GT HERE

	bcf_translate(out_hdr, hdr_in, line);
	if ( nswap[uphalf] ) phase_update(hdr_in, line, uphalf);
	//remove_info(out_hdr,line);
	//remove_format(out_hdr,line);
	if (bcf_write(fd, out_hdr, line) ) vrb.error("Failed to write the record output to file");
}

void concat::scan_overlap(const int ifname, const char * seek_chr, int seek_pos)
{
	bcf_srs_t * sr =  bcf_sr_init();
	sr->require_index = 1;
	sr->collapse = COLLAPSE_NONE;
	//sr->max_unpack = BCF_UN_FMT;

	int n_threads = options["threads"].as < int > ();
	if (n_threads > 1) bcf_sr_set_threads(sr, n_threads);

	if (!bcf_sr_add_reader (sr, filenames[ifname-2].c_str())) vrb.error("Problem opening/creating index file for [" + filenames[ifname-2] + "]");
	if (!bcf_sr_add_reader (sr, filenames[ifname-1].c_str())) vrb.error("Problem opening/creating index file for [" + filenames[ifname-1] + "]");

	int nset = 0;
	int n_sites_buff = 0;
	int n_sites_tot = 0;
	int last_pos=seek_pos;

	bcf1_t * line0 = NULL, * line1 = NULL;
	bcf_sr_seek(sr, seek_chr, seek_pos);
	while ((nset = bcf_sr_next_line (sr)))
	{
		if (nset==1)
		{
			if ( !bcf_sr_has_line(sr,0) && bcf_sr_region_done(sr,0)) break;  // no input from the first reader
			++n_sites_tot;
			continue;
		}
		line0 =  bcf_sr_get_line(sr, 0);
		if (line0->n_allele != 2) continue;

		line1 =  bcf_sr_get_line(sr, 1);

		//TODO: next GT HERE
		int nGTsa = bcf_get_genotypes(sr->readers[0].header, line0, &GTa, &mGTa);
		int nGTsb = bcf_get_genotypes(sr->readers[1].header, line1, &GTb, &mGTb);
		if ( nGTsa <= 0 || nGTsb <= 0 )
			vrb.error("GT field is not present in overlap at position: " + std::to_string(line0->pos + 1));

		update_distances();
		last_pos = line0->pos;
		++n_sites_buff;
		++n_sites_tot;
	}
	bcf_sr_destroy(sr);

	stats1D stats_all;
	stats1D phaseq;

	nswap[1]=0;
	for (int i = 0 ; i < nsamples; i++)
	{
		swap_phase[1][i] = nmatch[i] < nmism[i];
		nswap[1] += swap_phase[1][i];

		stats_all.push(nmatch[i] + nmism[i]);

		float f = 99;
        if ( nmatch[i] && nmism[i] )
        {
            // Entropy-inspired quality. The factor 0.7 shifts and scales to (0,1)
           float f0 = (float)nmatch[i]/(nmatch[i]+nmism[i]);
           f = (99*(0.7 + f0*logf(f0) + (1-f0)*logf(1-f0))/0.7);
        }
        phaseq.push(f);

		nmatch[i] = 0;
		nmism[i]  = 0;
	}

	if (n_sites_buff <=0) vrb.error("Overlap is empty");
	nsites_buff_d2.push_back(n_sites_buff/2);
	vrb.print("Buf " + stb.str(nsites_buff_d2.size() -1) + " ["+std::string(seek_chr)+":"+stb.str(seek_pos+1)+"-"+stb.str(last_pos+1)+"] [L_isec=" + stb.str(n_sites_buff) + " / L_tot=" + stb.str(n_sites_tot) + "] [Avg #hets=" + stb.str(stats_all.mean()) + "] [Switch rate=" + stb.str(nswap[1]*1.0 / nsamples) + "] [Avg phaseQ=" + stb.str(phaseq.mean()) + "]");
}
