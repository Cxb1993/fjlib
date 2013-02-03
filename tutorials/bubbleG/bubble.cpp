// ENABLE OR DISABLE PARALLEL SOLVER
// BY COMMENT LINE IN fjapp_bubbleNS.h"
// PARALLEL SOLVER ONLY WORKS IN LINUX
 
// ENABLE OR DISABLE LOGGING
//#define FJLIB_DISABLE_LOGGING

// DEFINE SOLVER HEADER HERE
#include "fjapp_bubbleNSG.h"

// DEFINE YOUR SOLVER HERE
typedef fjlib::TFJBubbleNSG_Solver my_Solver;

// DEFINED SURFACTANT #define _Surfactant;

#include "fjapp_bubble_driver.cpp"
