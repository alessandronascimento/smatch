################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/Coord.o \
../src/Gaussian.o \
../src/Mol.o \
../src/Optimization.o \
../src/Printer.o \
../src/SMatch.o 

CPP_SRCS += \
../src/Coord.cpp \
../src/Gaussian.cpp \
../src/Mol.cpp \
../src/Optimization.cpp \
../src/Printer.cpp \
../src/SMatch.cpp 

OBJS += \
./src/Coord.o \
./src/Gaussian.o \
./src/Mol.o \
./src/Optimization.o \
./src/Printer.o \
./src/SMatch.o 

CPP_DEPS += \
./src/Coord.d \
./src/Gaussian.d \
./src/Mol.d \
./src/Optimization.d \
./src/Printer.d \
./src/SMatch.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


