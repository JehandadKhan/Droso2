################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CKalmanFilter.cpp \
../src/CLasso.cpp \
../src/CMatrix.cpp \
../src/CMatrix3D.cpp \
../src/CRndStream.cpp \
../src/CSparseMatrix.cpp \
../src/CTimeVaryingNW.cpp \
../src/CVector.cpp \
../src/Kalman1.cpp \
../src/mmio.cpp 

OBJS += \
./src/CKalmanFilter.o \
./src/CLasso.o \
./src/CMatrix.o \
./src/CMatrix3D.o \
./src/CRndStream.o \
./src/CSparseMatrix.o \
./src/CTimeVaryingNW.o \
./src/CVector.o \
./src/Kalman1.o \
./src/mmio.o 

CPP_DEPS += \
./src/CKalmanFilter.d \
./src/CLasso.d \
./src/CMatrix.d \
./src/CMatrix3D.d \
./src/CRndStream.d \
./src/CSparseMatrix.d \
./src/CTimeVaryingNW.d \
./src/CVector.d \
./src/Kalman1.d \
./src/mmio.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	mpic++ -O3 -parallel -opt-matmul -mkl=parallel -openmp-report1 -par_report1 -vec-report1 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


