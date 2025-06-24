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

#ifndef _GTCHECK_ARGUMENT_SET_H
#define _GTCHECK_ARGUMENT_SET_H

#include "../utils/otools.h"
#include "../../versions/versions.h"

class gtcheck_argument_set {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

    static std::vector<std::string> mLegalOptions;

    uint32_t mSeed;
    uint32_t mNumThreads;
    std::vector<std::string> mInputFilenames;
    std::string mOutputFilename;

    bool mDeepCheck;

	gtcheck_argument_set(std::vector<std::string>& args)
    {
    	declare_options();
    	parse_command_line(args);
    	check_options();
    	verbose_files();
    	verbose_options();
    }

    void declare_options()
    {
    	try
    	{
        	bpo::options_description opt_base ("Basic options");
        	opt_base.add_options()
        			("help", "Produce help message")
					("threads", boost::program_options::value<uint32_t>(&mNumThreads)->default_value(1), "Number of threads.")
					("seed", boost::program_options::value<uint32_t>(&mSeed)->default_value(42), "Seed for RNG.")
					;

        	bpo::options_description opt_input("Input files");
        	opt_input.add_options()
        	    ("input,i", bpo::value<std::vector<std::string>>(&mInputFilenames)->multitoken()->required(),
        	     "Input genotype data in XCF format (e.g., --input file1.xcf file2.xcf)");

        	bpo::options_description opt_par("Parameters");
        	opt_input.add_options()
        	    ("deep-check", "Weather checking actual GTs instead of only GT counts. It does ignore the phase and missing GTs.");

        	bpo::options_description opt_output ("Output files");
        	opt_output.add_options()
        			("output,o", bpo::value< std::string >(&mOutputFilename), "Output file in XCF format")
        			("log", bpo::value< std::string >(), "Output log file")
					;

        	descriptions.add(opt_base).add(opt_input).add(opt_par).add(opt_output);
    	} catch ( const std::exception& e )
    	{
    		vrb.print("Available options:\n");
    		std::cout << descriptions << std::endl;
    		vrb.error(e.what());
    	}
    }

    void parse_command_line(std::vector < std::string > & args)
    {
    	try {
    		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
    		bpo::notify(options);
    	} catch ( const boost::program_options::error& e ) { std::cerr << "Error parsing command line arguments: " << std::string(e.what()) << std::endl; exit(0); }

    	if (options.count("l") && !vrb.open_log(options["l"].as < std::string > ()))
    		vrb.error("Impossible to create log file [" + options["l"].as < std::string > () +"]");

    	vrb.title("[XCFtools] Fill tags from/to XCF files");
    	vrb.bullet("Authors       : Olivier DELANEAU and Simone RUBINACCI");
    	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
    	vrb.bullet("Version       : 0." + std::string(XCFTLS_VERSION) + " / commit = " + std::string(__COMMIT_ID__) + " / release = " + std::string (__COMMIT_DATE__));
    	vrb.bullet("Run date      : " + tac.date());

    	if (options.count("help")) { std::cout << descriptions << std::endl; exit(0); }
    }

    void check_options()
    {
    	if (options.count("input") && options["input"].as < std::vector<std::string> > ().size() != 2)
			vrb.error("You must specify exactly two input XCF files with --input");

    	if (!options.count("output"))
    		vrb.error("You must specify an output XCF file with --output");

    	if (options.count("seed") && options["seed"].as < uint32_t > () < 0)
    		vrb.error("Random number generator needs a positive seed value");

    	if (options.count("threads") && options["threads"].as < uint32_t > () < 1)
    		vrb.error("You must use at least 1 thread");
    	mDeepCheck = options.count("deep-check");
    }

    void verbose_files() const
    {
    	std::array<std::string,2> no_yes = {"NO","YES"};

    	vrb.title("Files:");
    	vrb.bullet("Input files  : [" + mInputFilenames[0] + "] and [" + mInputFilenames[1] + "]");
    	vrb.bullet("Output       : [" + mOutputFilename + "]");
    	vrb.bullet("Deep check  : [" + no_yes[mDeepCheck] + "]");

    	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < std::string > () + "]");
    }

    void verbose_options() const
    {
    	std::array<std::string,2> no_yes = {"NO","YES"};

    	vrb.title("Parameters: ");
    	vrb.title("Other parameters");
    	vrb.bullet("Seed                : [" + stb.str(mSeed) + "]");
    	vrb.bullet("#Threads            : [" + stb.str(mNumThreads) + "]");
    }
};

#endif
