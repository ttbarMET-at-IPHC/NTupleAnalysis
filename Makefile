ObjSuf        = o
SrcSuf        = cc
ExeSuf        = 
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

ROOTCFLAGS    = $(shell root-config --cflags)
#ROOTLIBS      = $(shell root-config --libs) -lMinuit -lMathMore 
#ROOTGLIBS     = $(shell root-config --glibs) -lMinuit -lMathMore
ROOTLIBS      = -lGenVector -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lMinuit2 -lPhysics -lMathCore -lMathMore -lThread -pthread -lm -ldl -rdynamic -lTMVA
ROOTGLIBS     = -lGenVector -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMinuit -lMinuit2 -lPhysics -lMathCore -lMathMore -lThread -pthread -lm -ldl -rdynamic -lTMVA

# Linux with egcs
DEFINES       = -DNO_ORCA_CLASSES
CXX           = g++
CXXFLAGS	= -O -Wall -fPIC $(DEFINES)
LD		= g++
LDFLAGS		= -g -O -Wall -fPIC
SOFLAGS		= -shared

CXXFLAGS	+= $(ROOTCFLAGS) -I./ -I../../IPHCDataFormat/ -I../../
LIBS		= $(ROOTLIBS)  -lEG 
GLIBS		= $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard Selection/src/*.cc BckgdEstimation/src/*.cc Plots/src/*.cc Measurements/src/*.cc Tools/src/*.cc tinyxml/*.cc EffEstimation/src/*cc EventReco/src/*.cc JEC/src/*.cc BTagReshaping/src/*.cc)
HEADERS		= $(wildcard Selection/interface/*.h BckgdEstimation/interface/*.h Plots/interface/*.h Measurements/interface/*.h Tools/interface/*.h tinyxml/*.h EffEstimation/interface/*h EventReco/interface/*.h JEC/interface/*.h BTagReshaping/interface/*.cc)
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))


all:  libNTupleAna.so; cp  libNTupleAna.so .lib/libNTupleAna_`date +"%d-%m-%y_%H-%M-%S"`.so ; 

testDir:
	if [ -d ~/libs ]  then echo "" else mkdir ~/libs  fi

clean:
	@echo "Cleaning..."
	@rm -f $(OBJECTS) $(DEPENDS) *Dict.* core 
	make -C JetCorrections clean


.SUFFIXES: .$(SrcSuf) .C .o .so


libNTupleAna.so: $(OBJECTS) JetCorrections/libJetCorrections.so
	@echo "Building libNTupleAna..."
	$(LD) -L${ROOTSYS}/lib $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

JetCorrections/libJetCorrections.so:
	make -C JetCorrections lib
