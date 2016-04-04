TARGET = EP1

SRC = ./src/
OBJ = ./obj/

MGHSUB    = INITPT OBJFCN GRDFCN HESFCN
MGHFILES += $(addprefix $(SRC), $(MGHSUB:=.f))

MODNAMES  = buscalinear mgh
MODFILES  = $(MODNAMES:=.mod)
MODFILES := $(addprefix $(OBJ), $(MODFILES))

OBJFILES  = $(MODNAMES:=.o)
OBJFILES += $(MGHSUB:=.o)
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

$(OBJ)%.o $(OBJ)%.mod: $(SRC)%.f95
	$(FC) -J$(OBJ) -c $< -o $(OBJ)$*.o

$(OBJ)%.o: $(SRC)%.f
	$(FC) -c $< -o $@
