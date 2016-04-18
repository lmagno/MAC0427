TARGET = EP1

SRC = ./src/
OBJ = ./obj/
TEX = ./tex/
AUX = ./aux/

MGHSUB    = INITPT OBJFCN GRDFCN HESFCN
MGHFILES += $(addprefix $(SRC), $(MGHSUB:=.f))

MODNAMES  = methods mgh cholesky trisys utils ls stats paraboloid rosenbrock
MODFILES  = $(MODNAMES:=.mod)
MODFILES := $(addprefix $(OBJ), $(MODFILES))

OBJFILES  = $(MODNAMES:=.o)
OBJFILES += $(MGHSUB:=.o)
OBJFILES := $(addprefix $(OBJ), $(OBJFILES))
OBJFILES += $(OBJ)$(TARGET).o

FC = gfortran
FFLAGS = -Wall -g
LFLAGS = -lblas

.PHONY: directories clean tex

$(TARGET): directories $(OBJFILES)
	$(FC) -I$(OBJ) $(FFLAGS) $(OBJFILES) $(LFLAGS) -o $(TARGET)

$(OBJ)$(TARGET).o: $(TARGET).f08 $(MODFILES) $(SRC)runmethods.f08
	$(FC) -J$(OBJ) -I$(OBJ) $(FFLAGS) -c $< -o $@

$(OBJ)methods.o $(OBJ)methods.mod: $(SRC)methods.f08 $(OBJ)cholesky.mod $(OBJ)trisys.mod
	$(FC) -J$(OBJ) -c $< -o $(OBJ)methods.o

$(OBJ)cholesky.o $(OBJ)cholesky.mod: $(SRC)cholesky.f08 $(OBJ)utils.mod $(OBJ)trisys.mod
	$(FC) -J$(OBJ) $(FFLAGS) -c $< -o $(OBJ)cholesky.o

$(OBJ)%.o $(OBJ)%.mod: $(SRC)%.f08
	$(FC) -J$(OBJ) $(FFLAGS) -c $< -o $(OBJ)$*.o

$(OBJ)%.o: $(SRC)%.f
	$(FC) -c $< -o $@

directories:
	mkdir -p $(OBJ)

clean:
	rm -rf $(TARGET) $(OBJ) $(AUX)

tex:
	mkdir -p $(AUX); \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)$(TARGET).tex && \
	cp $(TEX)$(TARGET).bib $(AUX) && cd $(AUX) && bibtex $(TARGET) && cd - && \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)$(TARGET).tex && \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)$(TARGET).tex && \
	mv $(AUX)$(TARGET).pdf ./
