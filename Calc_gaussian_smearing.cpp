#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <string>
#include "Calc_gaussian_smearing.h"
using namespace std;
using std::string;

//************************************************************************
//double pi = 2.0*acos(0);
//float Pi = 2.0*acos(0);
//double DegToRad = 2.0*acos(0)/180.0;
//************************************************************************

double Gaussian(double X, double X0, double sigma){

double variable1 = (X - X0)/sigma;
double variable2 = sqrt(4.0*acos(0))*sigma;

return exp( -0.5*pow( variable1,2) ) / variable2;
}
//************************************************************************

void Calc_gaussian_smearing(double sigma, double range, double nominal_value, double mean_value, int smear_points, int progressive_point, double Results[2]){

double effective_value = mean_value - range*(1.0 - 2.0*progressive_point/(1.0*smear_points - 1.0));
double gaussian_weight = Gaussian(effective_value, mean_value, sigma)*2.0*range/(1.0*smear_points-1.0);

if(smear_points == 1){
	effective_value = nominal_value;
	gaussian_weight = 1.0;}

Results[0] = effective_value;
Results[1] = gaussian_weight;

return;
}
