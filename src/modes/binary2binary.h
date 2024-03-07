#ifndef BINARY2BINARY_H_
#define BINARY2BINARY_H_

#include <utils/otools.h>
#include <containers/bitvector.h>
#include <objects/sparse_genotype.h>
#include <utils/xcf.h>


#define CONV_BCF_BG	0
#define CONV_BCF_BH	1
#define CONV_BCF_SG	2
#define CONV_BCF_SH	3

class binary2binary {
public:
	//PARAM
	bitvector binary_bit_buf;
	std::vector<int32_t> sparse_int_buf;

	std::string region;
	int nthreads;
	int mode;
	float minmaf;
	bool drop_info;

	//CONSTRUCTORS/DESCTRUCTORS
	binary2binary(std::string, float, int, int, bool);
	virtual ~binary2binary();

	//PROCESS
	void convert(std::string, std::string);
	void convert(std::string, std::string, const bool exclude, const bool isforce, std::vector<std::string>& smpls);
	int32_t parse_genotypes(xcf_reader& XR, const uint32_t idx_file);


};

#endif /* BINARY2BINARY_H_ */
