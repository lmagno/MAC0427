SRC = ./src/
OBJ = ./obj/
TEX = ./tex/
AUX = ./aux/
DATA = ./data/

MGHSUB    = INITPT OBJFCN GRDFCN HESFCN
MGHFILES += $(addprefix $(SRC), $(MGHSUB:=.f))

MODNAMES  = methods mgh cholesky trisys utils ls stats paraboloid rosenbrock
MODFILES  = $(MODNAMES:=.mod)
MODFILES := $(addprefix $(OBJ), $(MODFILES))

OBJFILES  = $(MODNAMES:=.o)
OBJFILES += $(MGHSUB:=.o)
OBJFILES := $(addprefix $(OBJ), $(OBJFILES))
OBJFILES += $(OBJ)EP1.o

FC = gfortran
FFLAGS = -Wall -g
LFLAGS = -lblas

.PHONY: directories clean tex tex2

EP1: directories $(OBJFILES)
	$(FC) -I$(OBJ) $(FFLAGS) $(OBJFILES) $(LFLAGS) -o EP1

$(OBJ)EP1.o: $(SRC)EP1.f08 $(MODFILES) $(SRC)runmethods.f08
	$(FC) -J$(OBJ) -I$(OBJ) $(FFLAGS) -c $< -o $@

EP2: directories $(OBJ)EP2.o $(OBJ)constrained.o $(OBJ)test.o $(OBJ)stats2.o
	$(FC) -I$(OBJ) $(FFLAGS) $(OBJ)EP2.o $(OBJ)constrained.o $(OBJ)test.o $(OBJ)stats2.o $(OBJ)cholesky.o $(OBJ)trisys.o $(OBJ)ls.o $(LFLAGS) -o EP2

$(OBJ)EP2.o: $(SRC)EP2.f08 $(OBJ)constrained.mod $(OBJ)test.mod $(OBJ)stats2.mod
	$(FC) -J$(OBJ) -I$(OBJ) $(FFLAGS) -c $< -o $@

$(OBJ)constrained.o $(OBJ)constrained.mod: $(SRC)constrained.f08 $(OBJ)cholesky.mod $(OBJ)trisys.mod $(OBJ)ls.mod
	$(FC) -J$(OBJ) -c $< -o $(OBJ)constrained.o

$(OBJ)methods.o $(OBJ)methods.mod: $(SRC)methods.f08 $(OBJ)cholesky.mod $(OBJ)trisys.mod
	$(FC) -J$(OBJ) -c $< -o $(OBJ)methods.o

$(OBJ)cholesky.o $(OBJ)cholesky.mod: $(SRC)cholesky.f08 $(OBJ)utils.mod $(OBJ)trisys.mod
	$(FC) -J$(OBJ) $(FFLAGS) -c $< -o $(OBJ)cholesky.o

$(OBJ)%.o $(OBJ)%.mod: $(SRC)%.f08
	$(FC) -J$(OBJ) $(FFLAGS) -c $< -o $(OBJ)$*.o

$(OBJ)%.o: $(SRC)%.f
	$(FC) -c $< -o $@

tex:
	mkdir -p $(AUX); \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)EP1.tex && \
	cp $(TEX)EP1.bib $(AUX) && cd $(AUX) && bibtex EP1 && cd - && \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)EP1.tex && \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)EP1.tex && \
	mv $(AUX)EP1.pdf ./

TESTS := teste1 teste3 teste4
TESTS := $(addprefix $(AUX), $(TESTS))

tex2: $(TESTS:=.tex) $(TESTS:=.eps)
	cp $(TEX)EP2.bib $(AUX); \
	rubber -d --unsafe --into $(AUX) $(TEX)EP2.tex; \
	mv $(AUX)EP2.pdf .

$(AUX)%.tex $(AUX)%.eps: $(DATA)%.gp $(DATA)%.dat
	cd $(DATA); \
	gnuplot $*.gp; \
	mv $*.tex $*.eps .$(AUX)

directories:
	mkdir -p $(OBJ)

clean:
	rm -rf EP1 $(OBJ) $(AUX)*
