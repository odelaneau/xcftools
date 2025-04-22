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

#include <concat/concat_header.h>
#include <utils/xcf.h>

void concat::read_files_and_initialise() {
	//step0: Initialize seed & other
	rng.setSeed(options["seed"].as < int > ());

	//step1:read filenames
	std::string buffer;
	std::string filelist = options["input"].as < std::string > ();
	vrb.title("Read filenames in [" + filelist + "]");
	input_file fd(filelist);
	while (getline(fd, buffer))
	{
		filenames.push_back(buffer);
		if (!options.count("naive"))
		{
			xcf_reader XR(1);
			int32_t idx_file = XR.addFile(buffer);

			//Get file type
			int32_t type = XR.typeFile(idx_file);
			if (type != FILE_BINARY) vrb.error("[" + buffer + "] is not a XCF file");
			XR.close();
		}
	}
	vrb.bullet("#files = " + stb.str(filenames.size()));
	if (filenames.size() == 0) vrb.error("No filenames in input file.");

	nfiles = filenames.size();
}
