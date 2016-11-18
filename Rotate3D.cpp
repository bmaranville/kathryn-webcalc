#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <string>
#include "Rotate3D.h"
using namespace std;
using std::string;

//************************************************************************
//Example Usage:

//const int theta_steps = 45;
//double theta_step_size = 2.0;
//double theta_offset = 1.0;
//const int phi_steps = 9;
//double phi_step_size = 45.0;
//double phi_offset = 0.0;
//const int omega_steps = 7;//7 for 3Dobjects and 1 for linear objects
//double omega_step_size = 15.0;
//double omega_offset = 0.0;
//double RotationNorm[theta_steps][phi_steps][omega_steps];

//Optional:
//double XtoXTrans[theta_steps][phi_steps][omega_steps];
//double XtoYTrans[theta_steps][phi_steps][omega_steps];
//double XtoZTrans[theta_steps][phi_steps][omega_steps];
//double YtoXTrans[theta_steps][phi_steps][omega_steps];
//double YtoYTrans[theta_steps][phi_steps][omega_steps];
//double YtoZTrans[theta_steps][phi_steps][omega_steps];
//double ZtoXTrans[theta_steps][phi_steps][omega_steps];
//double ZtoYTrans[theta_steps][phi_steps][omega_steps];
//double ZtoZTrans[theta_steps][phi_steps][omega_steps];
//double CenterlineProjX[theta_steps][phi_steps];
//double CenterlineProjY[theta_steps][phi_steps];
//double CenterlineProjZ[theta_steps][phi_steps];

//for(int a=0; a<theta_steps; a++){
	//for(int b=0; b<phi_steps; b++){
	//double theta = (a*theta_step_size+theta_offset);
        //double phi = (b*phi_step_size+phi_offset);
		//for(int d=0; d<omega_steps; d++){//first d loop, omega loop
		//double omega = (d*omega_step_size + omega_offset);


		//Heart of sub-routine call:
		//double Trans_Matrix[3][3];
		//**with X=0, Y=1, and Z=2, wthis matrix computes  [start object coordinates along X, Y, Z] to [rotated object coordinates in X, Y, Z]
		//double CenterlineProjection[3];
		//**which computes projection onto X, Y, and Z for centerline of object (defined by theta and phi, and independent of omega)
		//RotationNorm[a][b][d] = Rotate3D(viewing_angle, theta, phi, omega, Trans_Matrix, CenterlineProjection);

		//Optional:
		//XtoXTrans[a][b][d] = Trans_Matrix[0][0];
		//XtoYTrans[a][b][d] = Trans_Matrix[0][1];
		//XtoZTrans[a][b][d] = Trans_Matrix[0][2];
		//YtoXTrans[a][b][d] = Trans_Matrix[1][0];
		//YtoYTrans[a][b][d] = Trans_Matrix[1][1];
		//YtoZTrans[a][b][d] = Trans_Matrix[1][2];
		//ZtoXTrans[a][b][d] = Trans_Matrix[2][0];
		//ZtoYTrans[a][b][d] = Trans_Matrix[2][1];
		//ZtoZTrans[a][b][d] = Trans_Matrix[2][2];
		//CenterlineProjX[a][b] = CenterlineProjection[0];
		//CenterlineProjY[a][b] = CenterlineProjection[1];
		//CenterlineProjZ[a][b] = CenterlineProjection[2];

		//}//end first d (omega) loop
        //}//end first b (phi) loop
//}//end first a (theta) loop


//************************************************************************
double Rotate3D(double viewing_angle, double theta, double phi, double omega, double Trans_Matrix[3][3], double CenterlineProjection[3]){
//Angles are given in degrees.
//theta is angle from x-axis, phi rotates within y-z plane with a contant theta, and omega is "roll" (for 3D objects) about
//centerline defined by theta and phi.

//Calculate the normalization factor and rotation matrices,
//where alpha = viewing_angle, beta is perpendicular to alpha and resides in x-y plane, and gamma = alpha x beta.
//When finsihed, x_pre-rotation = XtoXTrans[a][b][d], XtoYTrans[a][b][d], XtoZTrans[a][b][d], 
//y_pre-rotation = YtoXTrans[a][b][d], YtoYTrans[a][b][d], YtoZTrans[a][b][d], 
//z_pre-rotation = ZtoXTrans[a][b][d], ZtoYTrans[a][b][d], ZtoZTrans[a][b][d].
//phi is used to precess the pre-rotated x-axis about the viewing_angle, while omega is used to spin the 
//pre-rotated y and z axis of the structure within the rotated version beta-gamma plane.

double Norm = 0.0;
double DegToRad = acos(0.0)/90.0;
double XtoXTrans, XtoYTrans, XtoZTrans, YtoXTrans, YtoYTrans, YtoZTrans, ZtoXTrans, ZtoYTrans, ZtoZTrans;

double alpha[3], beta[3], gamma[3];
alpha[0] = cos(viewing_angle*DegToRad); alpha[1] = sin(viewing_angle*DegToRad); alpha[2] = 0;
beta[0] = cos((viewing_angle + 90.0)*DegToRad); beta[1] = sin((viewing_angle+90.0)*DegToRad); beta[2] = 0;
gamma[0] = 0.0; gamma[1] = 0.0; gamma[2] = 1.0;
int points = 3;
double LHS_Matrix[3][3];
LHS_Matrix[0][0] = alpha[0]; LHS_Matrix[0][1] = alpha[1]; LHS_Matrix[0][2] = alpha[2];
LHS_Matrix[1][0] = beta[0]; LHS_Matrix[1][1] = beta[1]; LHS_Matrix[1][2] = beta[2];
LHS_Matrix[2][0] = gamma[0]; LHS_Matrix[2][1] = gamma[1]; LHS_Matrix[2][2] = gamma[2];
double RHS_Matrix[3] = {cos(theta*DegToRad), sin(theta*DegToRad)*cos(phi*DegToRad), sin(theta*DegToRad)*sin(phi*DegToRad)};

if(alpha[0] <= 0.01 && alpha[0] >= -0.01){
LHS_Matrix[1][0] = alpha[0]; LHS_Matrix[1][1] = alpha[1]; LHS_Matrix[1][2] = alpha[2];
LHS_Matrix[0][0] = beta[0]; LHS_Matrix[0][1] = beta[1]; LHS_Matrix[0][2] = beta[2];
RHS_Matrix[0] = sin(theta*DegToRad)*cos(phi*DegToRad);
RHS_Matrix[1] = cos(theta*DegToRad);
}

//First row Gauss-Jordan ellimination:
for(int col = 0; col<points-1; col++){
	for(int row = col+1; row<points; row++){
		for(int p=1; p<points-col; p++){
		LHS_Matrix[row][col+p] = LHS_Matrix[row][col+p] -(LHS_Matrix[row][col] / LHS_Matrix[col][col])*LHS_Matrix[col][col+p];
		}//end p loop
		RHS_Matrix[row] = RHS_Matrix[row] -(LHS_Matrix[row][col] / LHS_Matrix[col][col])*RHS_Matrix[col];
		LHS_Matrix[row][col] = LHS_Matrix[row][col] -(LHS_Matrix[row][col] / LHS_Matrix[col][col])*LHS_Matrix[col][col];
	}//end row loop
}//end col loop

//Additional rows of Gauss-Jordan ellimination:
for(int col = 1; col<points; col++){
	for(int row = 0; row<col; row++){
		for(int p=1; p<points-col; p++){
		LHS_Matrix[row][col+p] = LHS_Matrix[row][col+p] -(LHS_Matrix[row][col] / LHS_Matrix[col][col])*LHS_Matrix[col][col+p];
		}//end p loop
	RHS_Matrix[row] = RHS_Matrix[row] -(LHS_Matrix[row][col] / LHS_Matrix[col][col])*RHS_Matrix[col];
	LHS_Matrix[row][col] = LHS_Matrix[row][col] -(LHS_Matrix[row][col] / LHS_Matrix[col][col])*LHS_Matrix[col][col];
	}//end row loop
}//end col loop

double Alignment[3];
for(int p=0; p<points; p++){
Alignment[p] = RHS_Matrix[p] / LHS_Matrix[p][p];
}//end p loop

double Unit = sqrt(Alignment[0]*Alignment[0] + Alignment[1]*Alignment[1] + Alignment[2]*Alignment[2]);
if(Unit > 1.001 || Unit < 0.999){cout << "Unit is " << Unit << endl;}
CenterlineProjection[0] = Alignment[0]/Unit;
CenterlineProjection[1] = Alignment[1]/Unit;
CenterlineProjection[2] = Alignment[2]/Unit;
//********************************************

XtoXTrans = Alignment[0]/Unit;
XtoYTrans = Alignment[1]/Unit;
XtoZTrans = Alignment[2]/Unit;

//Need to define a vector perpendicular to CenterlineProjectin[3] with length 1.0:
Alignment[0] = -XtoYTrans; Alignment[1] = XtoXTrans; Alignment[2] = 0.0;
Unit = sqrt(Alignment[0]*Alignment[0] + Alignment[1]*Alignment[1] + Alignment[2]*Alignment[2]);
//Here we don't expect Unit to be unity, so it needs to be calculated for the purpose of normalization:
YtoXTrans = Alignment[0]/Unit;
YtoYTrans = Alignment[1]/Unit;
YtoZTrans = Alignment[2]/Unit;

double dot = YtoXTrans*XtoXTrans + YtoYTrans*XtoYTrans + YtoZTrans*XtoZTrans;
if( pow(dot,2) > 0.000001){cout << "Error0" << endl;}

//Define third (perpendicular vector) using (XtoXTrans, XtoYTrans, XtoZTrans) cross
//(YtoXTrans, YtoYTrans, YtoZTrans) => (ZtoXTrans, ZtoYTrans, ZtoZTrans):
Alignment[0] = XtoYTrans*YtoZTrans - XtoZTrans*YtoYTrans;
Alignment[1] = -(XtoXTrans*YtoZTrans - XtoZTrans*YtoXTrans);
Alignment[2] = XtoXTrans*YtoYTrans - XtoYTrans*YtoXTrans;
Unit = sqrt(Alignment[0]*Alignment[0] + Alignment[1]*Alignment[1] + Alignment[2]*Alignment[2]);
if(Unit > 1.001 || Unit < 0.999){cout << "Unit3 is " << Unit << endl;}
ZtoXTrans = Alignment[0]/Unit;
ZtoYTrans = Alignment[1]/Unit;
ZtoZTrans = Alignment[2]/Unit;

dot = YtoXTrans*ZtoXTrans + YtoYTrans*ZtoYTrans + YtoZTrans*ZtoZTrans;
if( pow(dot,2) > 0.000001){cout << "Error1" << endl;}
dot = XtoXTrans*ZtoXTrans + XtoYTrans*ZtoYTrans + XtoZTrans*ZtoZTrans;
if( pow(dot,2) > 0.000001){cout << "Error2" << endl;}

//Rotate YTrans and ZTrans by omega:
double Intermed_YTrans[3];
double Intermed_ZTrans[3];

Intermed_YTrans[0] = YtoXTrans*cos(omega*DegToRad) - ZtoXTrans*sin(omega*DegToRad);
Intermed_YTrans[1] = YtoYTrans*cos(omega*DegToRad) - ZtoYTrans*sin(omega*DegToRad);
Intermed_YTrans[2] = YtoZTrans*cos(omega*DegToRad) - ZtoZTrans*sin(omega*DegToRad);

Intermed_ZTrans[0] = ZtoXTrans*cos(omega*DegToRad) + YtoXTrans*sin(omega*DegToRad);
Intermed_ZTrans[1] = ZtoYTrans*cos(omega*DegToRad) + YtoYTrans*sin(omega*DegToRad);
Intermed_ZTrans[2] = ZtoZTrans*cos(omega*DegToRad) + YtoZTrans*sin(omega*DegToRad);

Unit = sqrt(Intermed_YTrans[0]*Intermed_YTrans[0] + Intermed_YTrans[1]*Intermed_YTrans[1] + Intermed_YTrans[2]*Intermed_YTrans[2]);
if(Unit > 1.001 || Unit < 0.999){cout << "Unit4 is " << Unit << endl;}
YtoXTrans = Intermed_YTrans[0]/Unit;
YtoYTrans = Intermed_YTrans[1]/Unit;
YtoZTrans = Intermed_YTrans[2]/Unit;

Unit = sqrt(Intermed_ZTrans[0]*Intermed_ZTrans[0] + Intermed_ZTrans[1]*Intermed_ZTrans[1] + Intermed_ZTrans[2]*Intermed_ZTrans[2]);
if(Unit > 1.001 || Unit < 0.999){cout << "Unit5 is " << Unit << endl;}
ZtoXTrans = Intermed_ZTrans[0]/Unit;
ZtoYTrans = Intermed_ZTrans[1]/Unit;
ZtoZTrans = Intermed_ZTrans[2]/Unit;

Trans_Matrix[0][0] = XtoXTrans; Trans_Matrix[0][1] = XtoYTrans; Trans_Matrix[0][2] = XtoZTrans;
Trans_Matrix[1][0] = YtoXTrans; Trans_Matrix[1][1] = YtoYTrans; Trans_Matrix[1][2] = YtoZTrans;
Trans_Matrix[2][0] = ZtoXTrans; Trans_Matrix[2][1] = ZtoYTrans; Trans_Matrix[2][2] = ZtoZTrans;
Norm += sqrt(pow(sin(theta*DegToRad),2));

return Norm;
}
//************************************************************************
