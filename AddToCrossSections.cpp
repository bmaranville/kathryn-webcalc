#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <string>
#include "AddToCrossSections.h"
using namespace std;
using std::string;

//************************************************************************
void AddToCrossSections(double viewing_angle, double amp_nucl_real, double amp_nucl_img, double amp_mag_a_real, double amp_mag_a_img, double amp_mag_b_real, double amp_mag_b_img, double amp_mag_c_real, double amp_mag_c_img, double anglewt, double Q_gaussian_wt, double Norm_Volume, double CrossSection[4]){

//a is for magnetic projection || H (usually X in SANS Titan magnet), with b x c = a (b is usually magnetic projection along Y and c is usually magnetic projection along Z in SANS Titan magnet)

double GammaMatrix[3][3];
double QA_hat= cos(viewing_angle*acos(0.0)/90.0);
double QB_hat= sin(viewing_angle*acos(0.0)/90.0);
double QC_hat= 0.0;
GammaMatrix[0][0] = 1.0 - QA_hat*QA_hat; GammaMatrix[0][1] = -QA_hat*QB_hat; GammaMatrix[0][2] = -QA_hat*QC_hat;
GammaMatrix[1][0] = 1.0 - QB_hat*QB_hat; GammaMatrix[1][1] = -QA_hat*QB_hat; GammaMatrix[1][2] = -QB_hat*QC_hat;
GammaMatrix[2][0] = 1.0 - QC_hat*QC_hat; GammaMatrix[2][1] = -QA_hat*QC_hat; GammaMatrix[2][2] = -QB_hat*QC_hat;

double A_real = GammaMatrix[0][0]*amp_mag_a_real + GammaMatrix[0][1]*amp_mag_b_real + GammaMatrix[0][2]*amp_mag_c_real;
double A_img = GammaMatrix[0][0]*amp_mag_a_img + GammaMatrix[0][1]*amp_mag_b_img + GammaMatrix[0][2]*amp_mag_c_img;
double B_real = GammaMatrix[1][0]*amp_mag_b_real + GammaMatrix[1][1]*amp_mag_a_real + GammaMatrix[1][2]*amp_mag_c_real;
double B_img = GammaMatrix[1][0]*amp_mag_b_img + GammaMatrix[1][1]*amp_mag_a_img + GammaMatrix[1][2]*amp_mag_c_img;
double C_real = GammaMatrix[2][0]*amp_mag_c_real + GammaMatrix[2][1]*amp_mag_a_real + GammaMatrix[2][2]*amp_mag_b_real;
double C_img = GammaMatrix[2][0]*amp_mag_c_img + GammaMatrix[2][1]*amp_mag_a_img + GammaMatrix[2][2]*amp_mag_b_img;

CrossSection[0] += 0.5*(1E8)*anglewt*Q_gaussian_wt*(pow(amp_nucl_real-A_real,2)+pow(amp_nucl_img-A_img,2))/(Norm_Volume);
CrossSection[1] += 0.5*(1E8)*anglewt*Q_gaussian_wt*(pow(-B_real+C_img,2)+pow(-B_img-C_real,2))/(Norm_Volume);
CrossSection[2] += 0.5*(1E8)*anglewt*Q_gaussian_wt*(pow(amp_nucl_real+A_real,2)+pow(amp_nucl_img+A_img,2))/(Norm_Volume);
CrossSection[3] += 0.5*(1E8)*anglewt*Q_gaussian_wt*(pow(-B_real-C_img,2)+pow(-B_img+C_real,2))/(Norm_Volume);

return;
}
//************************************************************************
