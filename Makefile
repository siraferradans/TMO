#Compiling options
CPP=g++ 

CPPFLAGS=-g -O3
incdir1=./inc
incdir2=/opt/local/include/ImageMagick
incdirFFTW=./fftw-2.1.5/fftw
#Linking options
INCLUDES=-I$(incdir2) -I$(incdir1) -I$(incdirFFTW) 

libdirFFTW = ./fftw-2.1.5/fftw
LIBS=-lMagickCore -lMagick++ -llcms2 -ltiff -lfreetype -ljpeg -L/opt/local/lib -lfontconfig -lXext -lSM -lICE -lX11 -lXt -L/opt/local/lib -llzma -lbz2 -lz -lm -Wl,-framework,OpenCL -lm -lpthread -lltdl -L$(libdirFFTW) -lrfftw -lfftw

#Directories
SRC= ./src
OBJ = ./obj
bindir=./

TARGETS=$(CM)
OBJECTS=$(OBJCM)

#Items to build
CM=cm
OBJCM=mainCM.o Imagen.o inout.o utiles.o gradientes.o fft.o

all: $(CM)

cm: $(addprefix $(OBJ)/,$(OBJCM)) 
	@echo "[l] Linking file $@"
	$(CPP) $^ -o $@ $(INCLUDES) $(LIBS)

$(OBJ)/%.o :  $(addprefix $(SRC)/,%.cpp)
	$(CPP) -o $@ $(INCLUDES) $(CPPFLAGS) -c $^


all: $(cm)

clean: 
	@echo "[e] Cleaning..."
	-rm $(TARGETS)
	-rm $(OBJ)/*.o

cleanall: clean


