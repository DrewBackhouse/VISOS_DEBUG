#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <malloc.h>
#include "NeutrinoPropagator.h"
#include "BargerPropagator.h"
#include "EarthDensity.h"
#include "mosc.h"
#include "mosc3.h"
using namespace std;

//Defining functions

void Propagate(double, double, double, double, double, double, double, double,int, double, double);


int main (void) {

  printf("DEBUG: Start of main function\n");

  //Defining variables
  
  double num0 = 0.297;                  // Using 'defult' values from Visos
  double num1 = 0.425;                  // "
  double num2 = 0.0215;                 // " 
  double num3 = 0.0000737;              // "
  double num4 = 0.00249;                // "
  double num5 = 248;                    // "
  double energy = 0.003;                // Using 'Juno' Values
  double totaldistance = 100;           // " (But change 53 to 100)
  int nutype = -1;                      // "
  double density = 2.6;                 // "
  double error = 0.1;                   // "

  //calling function

  Propagate(num0, num1, num2, num3, num4, num5, energy, totaldistance , nutype, density, error);

}

void Propagate(double num0, double num1, double num2, double num3, double num4, double num5, double energy, double totaldistance ,int nutype, double density, double error) 
{
  
  printf("DEBUG: Start of function\n");
  
  //vector<double> X;
  //vector<double> Y;
  //vector<double> D;
  const int npoint = 44000;
  double* B = (double*)malloc(npoint * 8 * 3);
  double* A = (double*)malloc(npoint * 8 * 3);

  //Test  BargerPropagator *bNunue = new BargerPropagator();
  //Test  bNunue->UseMassEigenstates(false);

  //set PMNS can't be global, otherwise will be rewritten by other instances of bNunue
  //Test const double non_doubled_angle = true;
  // Set up 6 oscpars
  //(num0)= sin^2(theta12)
  //(num1)= sin^2(theta23)
  //(num2)=sin^2(theta13)
  //(num3)= dm^2_21
  //(num4)= dm^2_32
  //(num5) = dcp
  //Test  bNunue->SetMNS( num0, num2, num1, num3, num4, num5, energy, non_doubled_angle, nutype);
  
  const double distPoint = totaldistance/44000;
  double distance = 0;
  
  for(int ii=0; ii<npoint && distance<=totaldistance; ii++){
    // Set PMNS and propagate with given distance
    //Test  bNunue->propagateLinear(nutype, distance, density);
    
    //probablity projection.
    const double t_pnue = 1.0/3.0; //test  bNunue->GetProb( nutype, 1);
      const double t_pnumu = 1.0/3.0; //Test bNunue->GetProb( nutype, 2);
      const double t_pnutau = 1.0/3.0; //Test bNunue->GetProb( nutype, 3);
      //Test non_doubled_angleconst double t_sum = t_pnue + t_pnumu + t_pnutau;

    const double radangle = -135 * (M_PI/180);
    const double cosa = cos(radangle);
    const double sina = sin(radangle);

    const double tmpX =  cosa * t_pnumu - sina * t_pnutau;
    const double tmpY =  sina * t_pnumu + cosa * t_pnutau + 1/sqrt(2);
     const double tmpZ =  t_pnue;

    const double t_X = tmpX;

    const double radangle2 = -atan(sqrt(2));

    const double cosa2 = cos(radangle2);
    const double sina2 = sin(radangle2);
    
    const double t_Y =  cosa2 * tmpY - sina2 * tmpZ;
    //const double t_Z =  sina2 * tmpY + cosa2 * tmpZ;
    //X.push_back(t_X);
    //Y.push_back(t_Y);
    //D.push_back(distance);
    B[ii]=t_X;
    B[ii+44000]=t_Y;
    A[ii+88000]=distance;
    distance += distPoint;
  };
    
  //-----Function that calculates moving average of B and ouputs averaged  array A-----
  
    //Initialise variables
    double sum =0;
    //const int WindowSize = totaldistance * error / 100;
    const int WindowSize = 2000;
    
    //For loop over elements of array, excluding the last number of elements corisponding to WindowSize

    //Averaging over x components
    for (int i= 0; i<= (44000-WindowSize); i++){
      sum = 0; // Reset sum to zero

      //For loop over window
      for (int j = i; j < i + WindowSize; j++){
	sum += B[j]; //Adds value to the sum
      }
      A[i+(WindowSize/2)] = sum / (double)WindowSize; //Add averges to array. Shifts values by half window size to preserve symetries.

      //if statemet to check for first element and repeat this value over null entries.
      if (i==0){
	for (int k=0; k < (WindowSize/2); k++){
	  A[k]= sum / (double)WindowSize;
	}
      }
      //if statemet to check for last element and repeat this value over null entries
      if (i==(44000-WindowSize)){
	for (int k=(44000-WindowSize/2); k < (44000); k++){
	  A[k]= sum / (double)WindowSize;
	  cout << endl << i << " element of B " << B[i] << endl << endl; //prints element of array to terminal.
	}
      }
    }

    //Averaging over y components
    for (int i= 44000; i<= (88000-WindowSize); i++){
      sum = 0; // Reset sum to zero

      //For loop over window
      for (int j = i; j < i + WindowSize; j++){
	sum += B[j]; //Adds value to the sum
      }
      A[i+(WindowSize/2)] = sum / (double)WindowSize; //Add averges to array. Shifts values by half window size to preserve symetries

       //if statemet to check for first element and repeat this value over null entries.
      if (i==44000){
	for (int k=44000; k < (44000 + WindowSize/2); k++){
	  A[k]= sum / (double)WindowSize;
	}
      }
       //if statemet to check for last element and repeat this value over null entries.
      if (i==(88000-WindowSize)){
	for (int k=(88000-WindowSize/2); k < (88000); k++){
	  A[k]= sum / (double)WindowSize;
	}
      }
    }
}
