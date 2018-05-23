CFLAGS=-Wall -std=c++11
CFLAGS+=-g
#CFLAGS+=-O3

all:	demo0

test:	test_demo0

test_demo%:	demo%.tmp
	diff -s $< $(subst tmp,out,$<)

demo%.tmp:	demo%
	./$< < $<.in > $@

demo%:	demo%.cpp ebi.h ebi.cpp
	g++ $(CFLAGS) -o $@ $(patsubst %.h,,$^)

clean:
	rm -f demo? demo?.tmp

