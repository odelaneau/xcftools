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

#ifndef _OLIVIER_TOOLS_H //redefined to this name due to XCF stuff
#define _OLIVIER_TOOLS_H

//INCLUDE STANDARD TEMPLATE LIBRARY USEFULL STUFFS (STL)
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <bit>
#include <bitset>
#include <set>
#include <tuple>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <string>
#include <exception>
#include <cassert>
#include <limits>
#include <cstdint>

//INCLUDE BOOST USEFULL STUFFS (BOOST)
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/align/aligned_allocator.hpp>
#include <boost/math/distributions/chi_squared.hpp>

//INCLUDE HTS LIBRARY
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/regidx.h>
extern "C" {
	#include <htslib/vcf_sweep.h>
	#include <htslib/synced_bcf_reader.h>
	#include <htslib/vcf.h>
	#include <htslib/vcfutils.h>
}

//INCLUDES BASE STUFFS
#include <utils/compressed_io.h>
#include <utils/random_number.h>
#include <utils/basic_stats.h>
#include <utils/basic_algos.h>
#include <utils/string_utils.h>
#include <utils/timer.h>
#include <utils/verbose.h>

//#ifdef __RMATH_LIB__
//	#include <Rmath.h>
//#endif

//TYPEDEFS
template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

//CONSTANTS
#define RARE_VARIANT_FREQ	0.001f
#define HAP_NUMBER			8
#define MAX_AMB				22
#define MOD30BITS			0x40000000

//MACROS
#define DIV2(v)		((v)>>1)
#define MOD2(v)		((v)&1)
#define DIVU(v, d) 	(v+(d-1))/d
#define ROUND8(x)	((x+7)>>3)<<3

//NAMESPACE
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;

//MAKE SOME TOOLS FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
	random_number_generator rng;	//Random number generator
	string_utils stb;				//String manipulation
	basic_algos alg;				//Basic algorithms
	verbose vrb;					//Verbose
	timer tac;						//Timer
#else
	extern random_number_generator rng;
	extern string_utils stb;
	extern basic_algos alg;
	extern verbose vrb;
	extern timer tac;
#endif

#endif
