TARGET = EP1

SRC = ./
OBJ = ./obj/

MGHSRC    = ./MGH/
MGHSUB    = INITPT OBJFCN GRDFCN HESFCN
MGHFILES += $(addprefix $(MGHSRC), $(MGHSUB:=.f))

MODNAMES  = buscalinear mgh
MODFILES  = $(MODNAMES:=.mod)
OBJFILES  = $(MODNAMES:=.o)
OBJFILES += $(MGHSUB:=.o)
MODFILES := $(addprefix $(OBJ), $(MODFILES))
OBJFILES := $(addprefix $(OBJ), $(OBJFILES))
OBJFILES += $(OBJ)$(TARGET).o

FC = gfortran
FCLAGS = -Wall -Wimplicit-interface -g
LFLAGS = -lblas

.PHONY: directories clean

$(TARGET): directories $(OBJFILES)
	$(FC) -I$(OBJ) $(FCLAGS) $(OBJFILES) $(LFLAGS) -o $(TARGET)

directories:
	mkdir -p $(OBJ)

clean:
	rm -rf $(TARGET) $(OBJ)

$(OBJ)$(TARGET).o: $(TARGET).f95 $(MODFILES)
	$(FC) -I$(OBJ) -c $< -o $@

$(OBJ)mgh.mod $(OBJ)mgh.o: $(MGHSRC)mgh.f95
	$(FC) -J$(OBJ) -c $< -o $(OBJ)mgh.o

$(OBJ)buscalinear.mod $(OBJ)buscalinear.o: $(SRC)buscalinear.f95
	$(FC) -J$(OBJ) -c $< -o $(OBJ)buscalinear.o

$(OBJ)%.o: $(MGHSRC)%.f
	$(FC) -c $< -o $@
