/**
   \file
   Definition of SLT patterns and operations
   
   \author Luis Prado Jr (original code from M. Kleifges)
   \version $Id: SltPatternData.h 14717 2009-09-17 20:24:36Z lukas $
   \date 03 July 2004
*/

#ifndef _FdTriggerSimulatorOG_SltPatternData_h_
#define _FdTriggerSimulatorOG_SltPatternData_h_

static const char CVSId_FdTriggerSimulatorOG_SltPatternData[] = 
"$Id: SltPatternData.h 14717 2009-09-17 20:24:36Z lukas $";

namespace FdTriggerSimulatorOG {

  const short kSltRow[] = {-1,0,1,1,0,-1,-2,-2,-1,0,1,2,2,2,1,0,-1,-2};
  
  const short kSltColOdd[] = {1,1,1,0,-1,0,0,1,2,2,2,1,0,-1,-1,-2,-1,-1};
  
  const short kSltColEven[] = {0,1,0,-1,-1,-1,0,1,1,2,1,1,0,-1,-2,-2,-2,-1};
  
  const short kSltPattern[][3] = {{ 8, 1, 1},  { 1, 8, 1},  { 1, 1, 8},  
				  { 1, 1, 1},  {10, 2, 2},  { 2,10, 2},  
				  { 2, 2,10},  { 2, 2, 2},  {12, 3, 3},  
				  { 3,12, 3},  { 3, 3,12},  { 3, 3, 3},  
				  { 8, 1, 2},  { 1, 8, 2},  { 1, 1, 9},  
				  { 1, 1, 6},  {10, 2, 3},  { 2,10, 3},  
				  { 2, 2,11},  { 6, 8, 6},  {12, 3, 4},  
				  { 3,12, 4},  { 3, 3,13},  { 2, 3, 3},  
				  { 9, 1, 1},  { 2, 8, 1},  { 2, 1, 8},  
				  { 2, 1, 1},  {11, 2, 2},  { 3,10, 2},  
				  { 3, 2,10},  { 3, 2, 2},  {13, 3, 3},  
				  { 4,12, 3},  { 4, 3,12},  { 4, 3, 3},  
				  { 8, 1, 6},  { 1, 8, 6},  { 1, 1, 7},  
				  { 3, 4,13},  {10, 2, 1},  { 2,10, 1},  
				  { 2, 2, 9},  { 6, 1, 7},  {12, 3, 2},  
				  { 3,12, 2},  { 3, 3,11},  { 2, 3,11},  
				  { 7, 1, 1},  { 6, 8, 1},  { 6, 1, 8},  
				  { 6, 1, 1},  { 9, 2, 2},  { 1,10, 2},  
				  { 1, 2,10},  { 1, 2, 2},  {11, 3, 3},  
				  { 2,12, 3},  { 2, 3,12},  { 1,10, 1},  
				  { 8, 2, 1},  { 1, 9, 1},  { 2,12, 2},  
				  { 1, 1, 2},  {10, 3, 2},  { 2,11, 2},  
				  { 4, 3,13},  { 2, 2, 3},  {12, 4, 3},  
				  { 3,13, 3},  { 1, 2, 9},  { 3, 3, 4},  
				  { 9, 2, 1},  { 4,12, 4},  { 1, 2, 8},  
				  { 1, 2, 1},  { 3, 2,11},  { 3, 2, 3},  
				  { 2, 3,10},  { 2, 3, 2},  {13, 3, 4},  
				  { 3,11, 2},  { 3, 4,12},  { 3, 4, 3},  
				  { 8, 6, 1},  { 1, 7, 1},  { 2, 9, 1},  
				  {11, 3, 2},  {10, 1, 2},  { 2, 9, 2},  
				  {11, 2, 3},  { 2, 2, 1},  {12, 2, 3},  
				  { 3,11, 3},  { 2, 1, 9},  { 3, 3, 2},  
				  { 9, 1, 2},  { 3,10, 3},  { 1, 6, 8},  
				  { 1, 6, 1},  { 1, 9, 2},  { 2, 8, 2},  
				  { 2, 1,10},  { 2, 1, 2},  { 4,13, 3},  
				  { 2,11, 3},  { 3, 2,12},  { 3,13, 4} };
  
  const short kSltDiffOdd[][3] = {{20,41,40},  {21,41,40},  {21,20,40},  
				  {21,20,41},  {44,66,88},  {22,66,88},  
				  {22,44,88},  {22,44,66},  {24,47,48},  
				  {23,47,48},  {23,24,48},  {23,24,47},  
				  {20,41,63},  {21,41,63},  {21,20,63},  
				  {21,20,19},  {44,66,89},  {22,66,89},  
				  {22,44,89},  {23,20,19},  {24,47,26},  
				  {23,47,26},  {23,24,26},  {22,45,46},  
				  {43,42,63},  {22,42,63},  {22,43,63},  
				  {22,43,42},  {45,67,89},  {23,67,89},  
				  {23,45,89},  {23,45,67},  { 2,25,26},  
				  {21,24,47},  {21,23,47},  {21,23,24},  
				  {20,41,18},  {21,41,18},  {21,20,18},  
				  {23, 2, 4},  {44,66,87},  {22,66,87},  
				  {22,44,87},  {23,21,19},  {24,47,69},  
				  {23,47,69},  {23,24,69},  {22,45,68},  
				  {-2,19,18},  {23,20,41},  {23,21,41},  
				  {23,21,20},  {43,65,87},  {21,65,87},  
				  {21,43,87},  {21,43,65},  {45,46,69},  
				  {22,46,69},  {22,45,69},  {21,65,64},  
				  {20,42,63},  {21,42,63},  {22,46,68},  
				  {21,20,42},  {44,67,89},  {22,67,89},  
				  {21,23,25},  {22,44,67},  {24,25,26},  
				  {23,25,26},  {21,43,64},  {23,24,25},  
				  {43,65,64},  {21,24,25},  {21,43,63},  
				  {21,43,42},  {23,45,68},  {23,45,46},  
				  {22,45,89},  {22,45,67},  { 2,25, 4},  
				  {23,46,68},  {23, 2,26},  {23, 2,25},  
				  {20,19,18},  {21,19,18},  {22,65,64},  
				  {45,46,68},  {44,65,87},  {22,65,87},  
				  {45,67,68},  {22,44,65},  {24,46,69},  
				  {23,46,69},  {22,43,64},  {23,24,46},  
				  {43,42,64},  {23,67,68},  {21,-2,18},  
				  {21,-2,19},  {21,42,64},  {22,42,64},  
				  {22,43,87},  {22,43,65},  {21, 2,25},  
				  {22,67,68},  {23,45,69},  {23,25, 4} };
  
  const short kSltDiffEven[][3] = {{20,19,40},  {-1,19,40},  {-1,20,40},  
				   {-1,20,19},  {44,66,88},  {22,66,88},  
				   {22,44,88},  {22,44,66},  {24,25,48},  
				   { 1,25,48},  { 1,24,48},  { 1,24,25},  
				   {20,19,41},  {-1,19,41},  {-1,20,41},  
				   {-1,20,-3},  {44,66,67},  {22,66,67},  
				   {22,44,67},  { 1,20,-3},  {24,25,26},  
				   { 1,25,26},  { 1,24,26},  {22,23,46},  
				   {21,42,41},  {22,42,41},  {22,21,41},  
				   {22,21,42},  {23,45,67},  { 1,45,67},  
				   { 1,23,67},  { 1,23,45},  { 2, 3,26},  
				   {-1,24,25},  {-1, 1,25},  {-1, 1,24},  
				   {20,19,18},  {-1,19,18},  {-1,20,18},  
				   { 1, 2, 4},  {44,66,65},  {22,66,65},  
				   {22,44,65},  { 1,-1,-3},  {24,25,47},  
				   { 1,25,47},  { 1,24,47},  {22,23,68},  
				   {-2,-3,18},  { 1,20,19},  { 1,-1,19},  
				   { 1,-1,20},  {21,43,65},  {-1,43,65},  
				   {-1,21,65},  {-1,21,43},  {23,46,47},  
				   {22,46,47},  {22,23,47},  {-1,43,64},  
				   {20,42,41},  {-1,42,41},  {22,46,68},  
				   {-1,20,42},  {44,45,67},  {22,45,67},  
				   {-1, 1, 3},  {22,44,45},  {24, 3,26},  
				   { 1, 3,26},  {-1,21,64},  { 1,24, 3},  
				   {21,43,64},  {-1,24, 3},  {-1,21,41},  
				   {-1,21,42},  { 1,23,68},  { 1,23,46},  
				   {22,23,67},  {22,23,45},  { 2, 3, 4},  
				   { 1,46,68},  { 1, 2,26},  { 1, 2, 3},  
				   {20,-3,18},  {-1,-3,18},  {22,43,64},  
				   {23,46,68},  {44,43,65},  {22,43,65},  
				   {23,45,68},  {22,44,43},  {24,46,47},  
				   { 1,46,47},  {22,21,64},  { 1,24,46},  
				   {21,42,64},  { 1,45,68},  {-1,-2,18},  
				   {-1,-2,-3},  {-1,42,64},  {22,42,64},  
				   {22,21,65},  {22,21,43},  {-1, 2, 3},  
				   {22,45,68},  { 1,23,47},  { 1, 3, 4} };
  
  const short kSltMinRow[] = { 5, 5, 5, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4,
			       1, 1, 1, 4, 1, 1, 1, 1, 4, 4, 4, 3, 1, 1, 1, 1,
			       1, 2, 2, 2, 5, 5, 5, 1, 2, 2, 2, 4, 1, 1, 1, 1,
			       5, 4, 4, 3, 2, 2, 2, 2, 1, 1, 1, 3, 4, 4, 1, 3,
			       1, 1, 2, 1, 1, 1, 3, 1, 3, 2, 4, 3, 1, 1, 1, 1,
			       1, 1, 1, 1, 5, 5, 3, 1, 2, 2, 1, 2, 1, 1, 3, 1,
			       3, 1, 5, 4, 3, 3, 2, 2, 2, 1, 1, 1 };
  
  const short kSltMaxRow[] = {22,22,22,22,22,22,22,22,18,18,18,19,22,22,22,22,
			      21,21,21,21,18,18,18,20,22,22,22,22,21,21,21,21,
			      18,19,19,20,22,22,22,18,22,22,22,21,19,19,19,20,
			      22,21,21,21,22,22,22,22,19,19,19,22,22,22,20,22,
			      21,21,19,21,18,18,22,19,22,19,22,22,20,20,21,21,
			      18,20,18,19,22,22,22,20,22,22,20,22,19,19,22,20,
			      22,20,22,22,22,22,22,22,19,20,19,18 };
}

#endif  // _FdTriggerSimulatorOG_SltPatternData_h_
  

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "make -C .. -k"
// End: