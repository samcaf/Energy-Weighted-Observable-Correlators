#===========================================
# Install directory prefixes.
#===========================================
SOFTWARE_DIR=/Users/sam/Documents/Software

# Pythia
PYTHIA_BIN=$(SOFTWARE_DIR)/pythia8307/bin
PYTHIA_INCLUDE=$(SOFTWARE_DIR)/pythia8307/include
PYTHIA_LIB=$(SOFTWARE_DIR)/pythia8307/lib
PYTHIA_SHARE=$(SOFTWARE_DIR)/pythia8307/share/Pythia8
PYTHIA=$(PYTHIA_LIB)/libpythia8$(LIB_SUFFIX)

# Fastjet
FASTJET_LIB=$(SOFTWARE_DIR)/fastjet-3.4.0/src/.libs
FASTJET_INCLUDE=$(SOFTWARE_DIR)/fastjet-3.4.0/include
FASTJET=$(FASTJET_LIB)/libfastjet$(LIB_SUFFIX)


#===========================================
# Compilation flags
#===========================================
# Basic compilation
CXX=g++

CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_COMMON:=-I$(PYTHIA_INCLUDE) -I$(FASTJET_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PYTHIA_LIB) -L$(FASTJET_LIB) -Wl,-rpath,$(PYTHIA_LIB) -lpythia8 -lfastjet -ldl

CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib
