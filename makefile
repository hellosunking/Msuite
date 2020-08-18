Msuite: bin/preprocessor.pe bin/preprocessor.se bin/bowtie2.processer.pe bin/bowtie2.processer.se bin/rmdup.pe bin/rmdup.se bin/meth.caller.CpG bin/meth.caller.CpH bin/profile.DNAm.around.TSS util/bed2wig util/extract.meth.in.region
	@echo Build Msuite done.

cc=g++
## note that the g++ MUST support c++11 standard (i.e., version 4.8 or higher)
options=-std=c++11 -O2
multithread=-fopenmp #-pthread

bin/preprocessor.pe: src/preprocessor.pe.cpp src/common.h
	$(cc) $(options) $(multithread) -lz -o bin/preprocessor.pe src/preprocessor.pe.cpp

bin/preprocessor.se: src/preprocessor.se.cpp src/common.h
	$(cc) $(options) $(multithread) -o bin/preprocessor.se src/preprocessor.se.cpp

bin/bowtie2.processer.pe: src/bowtie2.processer.pe.cpp src/common.h
	$(cc) $(options) $(multithread) -o bin/bowtie2.processer.pe src/bowtie2.processer.pe.cpp

bin/bowtie2.processer.se: src/bowtie2.processer.se.cpp src/common.h
	$(cc) $(options) $(multithread) -o bin/bowtie2.processer.se src/bowtie2.processer.se.cpp

bin/rmdup.pe: src/rmdup.pe.cpp src/util.h src/util.cpp
	$(cc) $(options) -o bin/rmdup.pe src/rmdup.pe.cpp src/util.cpp

bin/rmdup.se: src/rmdup.se.cpp src/util.h src/util.cpp
	$(cc) $(options) -o bin/rmdup.se src/rmdup.se.cpp src/util.cpp

bin/meth.caller.CpG: src/meth.caller.CpG.cpp src/common.h src/util.h
	$(cc) $(options) -o bin/meth.caller.CpG src/meth.caller.CpG.cpp src/util.cpp

bin/meth.caller.CpH: src/meth.caller.CpH.cpp src/common.h src/util.h
	$(cc) $(options) -o bin/meth.caller.CpH src/meth.caller.CpH.cpp src/util.cpp

bin/profile.DNAm.around.TSS: src/profile.DNAm.around.TSS.cpp
	$(cc) $(options) -o bin/profile.DNAm.around.TSS src/profile.DNAm.around.TSS.cpp

util/bed2wig: util/bed2wig.cpp
	$(cc) $(options) -o util/bed2wig util/bed2wig.cpp

util/extract.meth.in.region: util/extract.meth.in.region.cpp
	$(cc) $(options) -o util/extract.meth.in.region util/extract.meth.in.region.cpp

clean:
	rm -f bin/preprocessor.pe bin/preprocessor.se bin/bowtie2.processer.pe bin/bowtie2.processer.se bin/rmdup.pe bin/rmdup.se bin/meth.caller util/bed2wig

