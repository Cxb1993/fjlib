// ENABLE OR DISABLE PARALLEL SOLVER
// BY COMMENT LINE IN fjapp_bubbleNS.h"
// PARALLEL SOLVER ONLY WORKS IN LINUX
 
// ENABLE OR DISABLE LOGGING
//#define FJLIB_DISABLE_LOGGING

// DEFINE SOLVER HEADER HERE
#include "fjapp_bubbleNSGD.h"

// DEFINE YOUR SOLVER HERE
typedef fjlib::TFJBubbleNSGD_Solver my_Solver;

// DEFINE SURFACTANT #define _Surfactant;
// DEFINE BULK_DIFFUSION #define _Bulk_Diffusion;

#include "fjapp_bubble_driver.cpp"
