CFLAGS=-Wall -std=c++11
CFLAGS+=-g
#CFLAGS+=-O3

all:	demo0.cpp ebi.h ebi.cpp
	g++ -std=c++11 -o demo demo0.cpp ebi.cpp

test:	test_demo0

test_demo%:	demo%.tmp
	@diff -s $< $(subst tmp,out,$<)

demo%.tmp:	demo%
	@if [ -f $<.in ] ; \
	then \
		./$< < $<.in > $@ ; \
	else \
		./$< > $@ ; \
	fi

demo%:	demo%.cpp ebi.h ebi.cpp
	g++ $(CFLAGS) -o $@ $(patsubst %.h,,$^)

clean:
	rm -f demo demo? demo?.tmp demo?? demo??.tmp

