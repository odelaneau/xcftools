#COMPILER MODE C++20
CXX=g++ -std=c++20


#create folders
dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

#COMPILER & LINKER FLAGS
CXXFLAGS+= -O3
LDFLAGS+= -O3
#CXXFLAGS+= -O0 -g
#LDFLAGS+= -O0

# Test if on x86 and target Haswell & newer.
# Disable this if building on x86 CPUs without AVX2 support.
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),x86_64)
    CXXFLAGS+= -march=x86-64-v3
endif

#COMMIT TRACING
COMMIT_VERS=$(shell git rev-parse --short HEAD)
COMMIT_DATE=$(shell git log -1 --format=%cd --date=short)
CXXFLAGS+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\"
CXXFLAGS+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"

# RMATH Support [YES/NO]
ifeq ($(RMATH_SUPPORT),)
	RMATH_SUPPORT=NO
endif
ifeq ($(RMATH_SUPPORT),YES)
	RMATH_INC=/usr/share/R/include/
	RMATH_LIB=/usr/lib/libRmath.a
#	DYN_LIBS+= -lzstd -lhts
	CXXFLAGS+= -D__RMATH_LIB__ -I$(RMATH_INC)
endif

# DYNAMIC LIBRARIES # Standard libraries are still dynamic in static exe
DYN_LIBS_FOR_STATIC=-lz -lpthread -lbz2 -llzma -lcrypto -ldeflate
# Non static exe links with all libraries
DYN_LIBS= -lboost_iostreams -lboost_program_options -lhts -pthread -lcurl -ldeflate

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

NAME=$(shell basename "$(CURDIR)")
BFILE=bin/$(NAME)
MACFILE=bin/$(NAME)_mac
EXEFILE=bin/$(NAME)_static

# Only search for libraries if goals != clean
ifeq (,$(filter clean,$(MAKECMDGOALS)))

#################################
# HTSLIB for static compilation #
#################################
# These are the default paths when installing htslib from source
HTSSRC=/usr/local
HTSLIB_INC=$(HTSSRC)/include/htslib
HTSLIB_LIB=$(HTSSRC)/lib/libhts.a

##########################################
# Boost libraries for static compilation #
##########################################
BOOST_INC=/usr/include

# If not set by user command, search for it
BOOST_LIB_IO?=$(shell whereis libboost_iostreams | grep -o '\S*\.a\b')
ifneq ($(suffix $(BOOST_LIB_IO)),.a)
    # If not found check default path
    ifeq ($(wildcard /usr/local/lib/libboost_iostreams.a),)
        # File does not exist
        $(warning libboost_iostreams.a not found, you can specify it with "make BOOST_LIB_IO=/path/to/lib...")
    else
        # File exists, set the variable
        BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
    endif
endif

# If not set by user command, search for it
BOOST_LIB_PO?=$(shell whereis libboost_program_options | grep -o '\S*\.a\b')
ifneq ($(suffix $(BOOST_LIB_PO)),.a)
    # If not found check default path
    ifeq ($(wildcard /usr/local/lib/libboost_program_options.a),)
        # File does not exist
        $(warning libboost_program_options.a not found, you can specify it with "make BOOST_LIB_PO=/path/to/lib...")
    else
        # File exists, set the variable
        BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
    endif
endif

# Endif makefile goals != clean
endif

#CONDITIONAL PATH DEFINITON
desktop: $(BFILE)

olivier: HTSSRC=$(HOME)/Tools
olivier: HTSLIB_INC=$(HTSSRC)/htslib-1.15
olivier: HTSLIB_LIB=$(HTSSRC)/htslib-1.15/libhts.a
olivier: BOOST_INC=/usr/include
olivier: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
olivier: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
olivier: $(BFILE)

rgc: HTSSRC=/mnt/efs_v2/agds_methods/users/olivier.delaneau/LIBS
rgc: HTSLIB_INC=$(HTSSRC)/htslib-1.18
rgc: HTSLIB_LIB=$(HTSSRC)/htslib-1.18/libhts.a
rgc: BOOST_INC=/usr/include
rgc: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
rgc: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
rgc: $(BFILE)

static_exe: HTSSRC=/home/srubinac/git
static_exe: HTSLIB_INC=$(HTSSRC)/htslib_minimal
static_exe: HTSLIB_LIB=$(HTSSRC)/htslib_minimal/libhts.a
static_exe: CXXFLAG=-O3 -mavx2 -mfma -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
static_exe: LDFLAG=-O3
static_exe: $(EXEFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAGS) $^ -o $@ $(DYN_LIBS)

$(EXEFILE): $(OFILE)
	$(CXX) $(LDFLAGS) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS_FOR_STATIC)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAGS) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE) $(EXEFILE)
