/*
 * Kalman1.h
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#ifndef KALMAN1_H_
#define KALMAN1_H_
#include <mkl.h>
#include <mkl_vsl.h>
#include <stdio.h>
#include <iostream>
//#include <math.h>
#include <cstdlib>
#include <assert.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <stdlib.h>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

#define eprintf(MSG) \
{													\
	int r = 0;									\
	MPI_Comm_rank(MPI_COMM_WORLD, &r);					\
	printf("[%d]%s: %d %s\n",r,__FILE__,__LINE__,MSG); \
}


//#include "log4cpp/Category.hh"
//#include "log4cpp/Configurator.hh"

//#include "CVector.h"
//#include "CMatrix.h"
//#include "CMatrix3D.h"




#endif /* KALMAN1_H_ */
