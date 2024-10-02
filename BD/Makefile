# имя программы:
TARGET=../BrownianDynamic
TARGET_CUBIC=../BrownianDynamic_Cubic
TARGET_CYLINDER=../BrownianDynamic_Cylinder
TARGET_HI=../BrownianDynamic_HI

#	BoundaryCondition.90 \
# исходники для сборки: 
SOURCES = \
	CommonModule.f90 TestErrors.f90 \
	BoxMuller.f90 \
	Surface.f90 Potential.f90 \
	PeriodicCondition.f90 \
	InputData.f90 OutputData.f90 \
	HydroDynamic.f90 AntiGyration.f90 \
	ForceEnergy.f90 \
	BrownianDynamic.f90 \

OBJECTS=$(SOURCES:%.f90=%.o)

TYPE=master

# про флаги:
# https://www.opennet.ru/docs/RUS/cpp/cpp-10.html
#### вместо раскоментированных #define в коде:
# для граничных условий:
CPPFLAGS_CYLINDER=BOUNDARY_CONDITION_CYLINDRICAL
CPPFLAGS_CUBIC=BOUNDARY_CONDITION_CUBIC
# для гидродинамических взаимодействий:
CPPFLAGS_HI=HYDRODYNAMIC_INTERACTION

# простая сборка:
all: $(TARGET)

$(OBJECTS): $(SOURCES)

$(TARGET): $(OBJECTS)
	gfortran -O3 -o $(TARGET) -cpp $(SOURCES) 
clean:
	$(RM) $(TARGET) *.mod

######### компиляция с граничными условиями:
# для кубической системы координат:
cubic: $(TARGET_CUBIC)

$(OBJECTS): $(SOURCES)

$(TARGET_CUBIC): $(OBJECTS)
	gfortran -O3 -o $(TARGET_CUBIC) -cpp  -D$(CPPFLAGS_CUBIC)=definition $(SOURCES)
	
# для цилиндрической системы координат:
cylinder: $(TARGET_CYLINDER)

$(OBJECTS): $(SOURCES)

$(TARGET_CYLINDER): $(OBJECTS)
	gfortran -O3 -o $(TARGET_CYLINDER) -cpp  -D$(CPPFLAGS_CYLINDER)=definition $(SOURCES) 

######### компиляция с гидродинамическими взаимодействиями:
hi: $(TARGET_HI)

$(OBJECTS): $(SOURCES)

$(TARGET_HI): $(OBJECTS)
	gfortran -O3 -o $(TARGET_HI) -cpp  -D$(CPPFLAGS_HI)=definition $(SOURCES) 

	
.PHONY: all clean