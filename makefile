#COMPILER MODE C++17
CXX=g++ -std=c++17


#create folders
dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

#COMPILER & LINKER FLAGS
CXXFLAG=-O3
LDFLAG=-O3

#CXXFLAG=-O0 -g -Wno-ignored-attributes
#LDFLAG=-O0

#COMMIT TRACING
COMMIT_VERS=$(shell git rev-parse --short HEAD)
COMMIT_DATE=$(shell git log -1 --format=%cd --date=short)
CXXFLAG+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\"
CXXFLAG+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"

# RMATH Support [YES/NO]
ifeq ($(RMATH_SUPPORT),)
	RMATH_SUPPORT=NO
endif
ifeq ($(RMATH_SUPPORT),YES)
	RMATH_INC=/usr/share/R/include/
	RMATH_LIB=/usr/lib/libRmath.a
#	DYN_LIBS+= -lzstd -lhts
	CXXFLAG+= -D__RMATH_LIB__ -I$(RMATH_INC)
endif

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lpthread -lbz2 -llzma -lcurl -lcrypto -ldeflate

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

NAME=$(shell basename "$(CURDIR)")
BFILE=bin/$(NAME)
MACFILE=bin/$(NAME)_mac
EXEFILE=bin/$(NAME)_static

#CONDITIONAL PATH DEFINITON
desktop: HTSSRC=/home/srubinac/git
desktop: HTSLIB_INC=$(HTSSRC)/htslib-1.17
desktop: HTSLIB_LIB=$(HTSSRC)/htslib-1.17/libhts.a
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
desktop: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
desktop: $(BFILE)

olivier: HTSSRC=$(HOME)/Tools
olivier: HTSLIB_INC=$(HTSSRC)/htslib-1.15
olivier: HTSLIB_LIB=$(HTSSRC)/htslib-1.15/libhts.a
olivier: BOOST_INC=/usr/include
olivier: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
olivier: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
olivier: $(BFILE)

static_exe: CXXFLAG=-O2 -mavx2 -mfma -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
static_exe: LDFLAG=-O2
static_exe: HTSSRC=/home/srubinac/git
static_exe: HTSLIB_INC=$(HTSSRC)/htslib_minimal
static_exe: HTSLIB_LIB=$(HTSSRC)/htslib_minimal/libhts.a
static_exe: BOOST_INC=/usr/include
static_exe: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
static_exe: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
static_exe: $(EXEFILE)

# static desktop Robin
static_exe_robin_desktop: CXXFLAG=-O2 -mavx2 -mfma -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
static_exe_robin_desktop: LDFLAG=-O2
static_exe_robin_desktop: HTSSRC=/home/robin/Dropbox/LIB
static_exe_robin_desktop: HTSLIB_INC=$(HTSSRC)/htslib_minimal
static_exe_robin_desktop: HTSLIB_LIB=$(HTSSRC)/htslib_minimal/libhts.a
static_exe_robin_desktop: BOOST_INC=/usr/include
static_exe_robin_desktop: BOOST_LIB_IO=$(HTSSRC)/boost/lib/libboost_iostreams.a
static_exe_robin_desktop: BOOST_LIB_PO=$(HTSSRC)/boost/lib/libboost_program_options.a
static_exe_robin_desktop: $(EXEFILE)

mac_apple_silicon: CXXFLAG=-O3 -mcpu=apple-m1 -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
mac_apple_silicon: LDFLAG=-O3
mac_apple_silicon: HTSLIB_LIB=-lhts
mac_apple_silicon: BOOST_LIB_IO=-lboost_iostreams
mac_apple_silicon: BOOST_LIB_PO=-lboost_program_options
mac_apple_silicon: $(MACFILE)

mac_apple_silicon_static: CXXFLAG=-O3 -mcpu=apple-m1 -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
mac_apple_silicon_static: LDFLAG=-O3
mac_apple_silicon_static: HTSLIB_LIB=/opt/homebrew/opt/htslib/lib/libhts.a
mac_apple_silicon_static: BOOST_LIB_IO=/opt/homebrew/opt/boost/lib/libboost_iostreams.a
mac_apple_silicon_static: BOOST_LIB_PO=/opt/homebrew/opt/boost/lib/libboost_program_options.a
mac_apple_silicon_static: $(MACFILE)

#COMPILATION RULES
all: desktop

$(MACFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

$(EXEFILE): $(OFILE)
	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE) $(EXEFILE)
