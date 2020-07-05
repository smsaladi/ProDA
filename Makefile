################################################################################
# Makefile for proda
################################################################################

CXX = g++

OTHERFLAGS = -DVERSION="\"1.00\""

CXXFLAGS = -O3 -W -Wall -pedantic -DNDEBUG $(OTHERFLAGS) -mmmx -msse -msse2 -mfpmath=sse -funroll-loops -fomit-frame-pointer 

################################################################################
# 3) Dependencies
################################################################################

TARGETS = proda
OBJECTS = AlignedFragment.o Assert.o Block.o Consistency.o GlobalAlign.o LocalAlign.o Main.o PairAligner.o Matrix.o MultiSequence.o ProbModel.o Score.o ScoreMatrix.o Sequence.o SparseMatrix.o Tree.o Utilities.o

.PHONY : all
all : $(TARGETS)

proda : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -lm $(OBJECTS) -o proda

Assert.o: Assert.h
AlignedFragment.o: AlignedFragment.h Utilities.h
Block.o: Block.h Sequence.h AlignedFragment.h Sequence.h Utilities.h Types.h
GlobalAlign.o: Assert.h GlobalAlign.h Matrix.h MultiSequence.h ProbModel.h Score.h Sequence.h
Main.o: Assert.h Block.h Consistency.h GlobalAlign.h Matrix.h MultiSequence.h ProbModel.h SparseMatrix.h Tree.h Utilities.h PairAligner.h
PairAligner.o: PairAligner.h Sequence.h AlignedFragment.h ProbModel.h
Matrix.o: Matrix.h Score.h ScoreMatrix.h SparseMatrix.h
MultiSequence.o: Assert.h MultiSequence.h Sequence.h Utilities.h
ProbModel.o: Matrix.h ProbModel.h Score.h
Score.o: Score.h
ScoreMatrix.o: ScoreMatrix.h Score.h
Sequence.o: Assert.h Sequence.h
SparseMatrix.o: Matrix.h Score.h SparseMatrix.h
Tree.o: GlobalAlign.h Matrix.h MultiSequence.h Tree.h 
Utilities.o: Assert.h Utilities.h

.PHONY : clean
clean:
	rm -f *.o
