PROGRAM       =  JetJt
version       = JTKT
CXX           = g++
#CXXFLAGS      = -O -Wall -g -Wno-deprecated -bind_at_load -D$(version)
CXXFLAGS      = -O -Wall -g -Wno-deprecated -fPIC -D$(version) #-ggdb
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)
CXXFLAGS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )
LDFLAGS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins ) 
INCS = -I/n/work00/towisnel/Herwig/Install/include -I/n/work00/towisnel/Herwig/test/src
#INCS  += -I$(HERWIG)/include -I$(THEPEG)/include
CXXFLAGS  += $(INCS)

HDRSDICT = src/AliJBaseTrack.h src/AliJCard.h src/AliJJet.h src/AliJBaseCard.h src/AliJJetJtHistos.h src/AliJHistManager.h src/AliJHistogramInterface.h
           
HDRS	+= $(HDRSDICT)  \
                        nanoDict.h


SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM) :     $(OBJS) src/AliJConst.h
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -fPIC -shared -lEG -lPhysics -L$(PWD) $(PROGRAM).cxx $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM).so
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) -fPIC $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl

cl:  clean $(PROGRAM)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -fPIC -c -D$(version) $(HDRSDICT)
