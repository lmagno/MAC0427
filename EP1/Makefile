TARGET = EP1

SRC = ./src/
OBJ = ./obj/
TEX = ./tex/
AUX = ./tex/aux/

MGHSUB    = INITPT OBJFCN GRDFCN HESFCN
MGHFILES += $(addprefix $(SRC), $(MGHSUB:=.f))

MODNAMES  = buscalinear mgh cholesky trisys utils ls stats
MODFILES  = $(MODNAMES:=.mod)
MODFILES := $(addprefix $(OBJ), $(MODFILES))

OBJFILES  = $(MODNAMES:=.o)
OBJFILES += $(MGHSUB:=.o)
OBJFILES := $(addprefix $(OBJ), $(OBJFILES))
OBJFILES += $(OBJ)$(TARGET).o

FC = gfortran
FCLAGS = -Wall -Wimplicit-interface -g
LFLAGS = -lblas

.PHONY: directories clean tex

$(TARGET): directories $(OBJFILES)
	$(FC) -I$(OBJ) $(FCLAGS) $(OBJFILES) $(LFLAGS) -o $(TARGET)

$(OBJ)$(TARGET).o: $(TARGET).f95 $(MODFILES)
	$(FC) -J$(OBJ) -I$(OBJ) -c $< -o $@

$(OBJ)buscalinear.o $(OBJ)buscalinear.mod: $(SRC)buscalinear.f95 $(OBJ)cholesky.mod $(OBJ)trisys.mod
	$(FC) -J$(OBJ) -c $< -o $(OBJ)buscalinear.o

$(OBJ)cholesky.o $(OBJ)cholesky.mod: $(SRC)cholesky.f95 $(OBJ)utils.mod $(OBJ)trisys.mod
	$(FC) -J$(OBJ) -c $< -o $(OBJ)cholesky.o

$(OBJ)%.o $(OBJ)%.mod: $(SRC)%.f95
	$(FC) -J$(OBJ) -c $< -o $(OBJ)$*.o

$(OBJ)%.o: $(SRC)%.f
	$(FC) -c $< -o $@

directories:
	mkdir -p $(OBJ)

clean:
	rm -rf $(TARGET) $(OBJ) $(AUX)

tex:
	mkdir -p $(AUX); \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)$(TARGET).tex && \
	pdflatex -halt-on-error -output-directory $(AUX) $(TEX)$(TARGET).tex && \
	mv $(AUX)$(TARGET).pdf ./
