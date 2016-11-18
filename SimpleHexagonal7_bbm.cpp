#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cstdio>
#include "Calc_gaussian_smearing.h"
#include "Calc_gaussian_smearing.cpp"
#include "Rotate3D.h"
#include "Rotate3D.cpp"
#include "AddToCrossSections.h"
#include "AddToCrossSections.cpp"

#include <emscripten/bind.h>
#include "json.hpp"

using namespace emscripten;
using json = nlohmann::json;

void fill_double_from_json(double* to_fill, json j) {
    //json j = json::parse(vector_json);
    int index = 0;
    // iterate the array
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
	  //cout << ((double) *it);
      to_fill[index++] = ((double) *it);
    }
    //cout << endl;
}

json default_params = json::array({
	{"correlation", 1.0}, // 1.0 = picks axis closest to X; 0 = picks random axis
    {"YesNo_q_smearing", 0}, //0 = no; 1=yes
    {"scatt_centers", 20}, //20, up to 90 possible
    {"structure_tilt_angle", 5.0}, //in x-y plane
    {"min_parallel_projection", 0.3}, 
    {"pinned_frac", 0.025},

    {"magcore_sld", 1.4E-6},
    {"magnetic_radius", 35.0},

    {"core_sld", 6.9E-6}, //6.9E-6 for Fe3O4
    {"shell1_sld", 5.2E-6}, //5.18E-6 for MnFe2O4
    {"shell2_sld", 1.5E-6}, //1.5E-6 for Mn3O4?
    {"solvent_sld", 0.0},
    {"core_radius", 10.0}, //10.0
    {"shell1_thickness", 8.0}, //18.0
    {"shell2_thickness", 20.0}, //8.0

    {"length", 84.0}, //85.0
    {"volume_fraction", 0.1}, //0.045 for 20 scatt centers
    {"viewing_angle", 0},

    {"background", 0.0}, //0.01,
    {"SFbackground", 0.05}, //0.05,

    {"theta_steps", 45}, //45.0 or 23
    {"theta_step_size", 2.0}, //2.0 4.0
    {"theta_offset", 1.0}, //1.0
    {"phi_steps", 45}, //45 or 30 or 18
    {"phi_step_size", 8.0}, //8.0 or 12.0 or 20.0
    {"phi_offset", 0.0}, //0.0
    {"omega_steps", 7}, //7
    {"omega_step_size", 15.0}, //15.0
    {"omega_offset", 0.0}
});

std::string get_default_params() {
	return default_params.dump();
};

//std::vector<std::string> get_default_params_keys() {
	
//*************************************************************************

//************************************************************************
std::string calculate(
	double correlation, //1.0 = picks axis closest to X; 0 = picks random axis
	double YesNo_q_smearing, //0 = no; 1=yes
	double scatt_centers, //20, up to 90 possible (convert to int)
	double structure_tilt_angle, //in x-y plane
	double min_parallel_projection, 
	double pinned_frac, 

	double magcore_sld,
	double magnetic_radius,

	double core_sld, //6.9E-6 for Fe3O4
	double shell1_sld, //5.18E-6 for MnFe2O4
	double shell2_sld, //1.5E-6 for Mn3O4?
	double solvent_sld,
	double core_radius, //10.0
	double shell1_thickness, //18.0
	double shell2_thickness, //8.0

	double length, //85.0
	double volume_fraction, //0.045 for 20 scatt centers
	double viewing_angle, 

	double background, //0.01;
	double SFbackground, //0.05;

	double theta_steps_d, //45.0 or 23
	double theta_step_size, //2.0 4.0
	double theta_offset, //1.0
	double phi_steps_d, //45 or 30 or 18
	double phi_step_size, //8.0 or 12.0 or 20.0
	double phi_offset, //0.0
	double omega_steps_d, //7
	double omega_step_size, //15.0
	double omega_offset,
	
	std::vector<double> Data_Q,
	std::vector<double> Data_DeltaQ,
	std::vector<double> Data_MeanQ
	
)
{
//*************************************************************************
//User Inputs:
//*************************************************************************
const int theta_steps = (int) theta_steps_d;
const int phi_steps = (int) phi_steps_d;
const int omega_steps = (int) omega_steps_d;
//*************************************************************************
//Non-user defined inputs:
//*************************************************************************
double PiNum = 2.0*acos(0.0);

//*************************************************************************
//Get SANS data for template and Q-smearing:
//*************************************************************************

//double Data_Q[1000], Data_I[1000], Data_DeltaI[1000], Data_DeltaQ[1000], Data_MeanQ[1000], Data_Shadow[1000];
//double Data_I[1000], Data_DeltaI[1000], Data_DeltaQ[1000], Data_MeanQ[1000], Data_Shadow[1000];
//double Data_Q2[1000], Data_I2[1000], Data_DeltaI2[1000], Data_DeltaQ2[1000], Data_MeanQ2[1000], Data_Shadow2[1000];
//double Data_Q3[1000], Data_I3[1000], Data_DeltaI3[1000], Data_DeltaQ3[1000], Data_MeanQ3[1000], Data_Shadow3[1000];
//double Data_Q4[1000], Data_I4[1000], Data_DeltaI4[1000], Data_DeltaQ4[1000], Data_MeanQ4[1000], Data_Shadow4[1000];

int remove_points = 25;//40

const int QPoints = Data_Q.size();
//std::copy(Q_values.begin(), Q_values.end(), Data_Q);

//Get_6ColumnSANSdata(remove_points, FileToReadIn, Data_Q, Data_I, Data_DeltaI, Data_DeltaQ, Data_MeanQ, Data_Shadow);
//Get_6ColumnSANSdata(remove_points, FileToReadIn2, Data_Q2, Data_I2, Data_DeltaI2, Data_DeltaQ2, Data_MeanQ2, Data_Shadow2);
//Get_6ColumnSANSdata(remove_points, FileToReadIn3, Data_Q3, Data_I3, Data_DeltaI3, Data_DeltaQ3, Data_MeanQ3, Data_Shadow3);
//Get_6ColumnSANSdata(remove_points, FileToReadIn4, Data_Q4, Data_I4, Data_DeltaI4, Data_DeltaQ4, Data_MeanQ4, Data_Shadow4);

double Model1[1000];
double Model2[1000];
double Model3[1000];
double Model4[1000];

//*************************************************************************
//Define Scattering Centers and Normalization Volume:
//*************************************************************************

double coordinates[90][3] = {0.0};
double Rcoordinates[90][3] = {0.0};
Rcoordinates[0][0] = 0.0; Rcoordinates[0][1] = 0.0;
Rcoordinates[1][0] = -length; Rcoordinates[1][1] = 0.0;
Rcoordinates[2][0] = length; Rcoordinates[2][1] = 0.0;
Rcoordinates[3][0] = length/2.0; Rcoordinates[3][1] = sqrt(3.0)*length/2.0;
Rcoordinates[4][0] = -length/2.0; Rcoordinates[4][1] = -sqrt(3.0)*length/2.0;
Rcoordinates[5][0] = -length/2.0; Rcoordinates[5][1] = sqrt(3.0)*length/2.0;
Rcoordinates[6][0] = length/2.0; Rcoordinates[6][1] = -sqrt(3.0)*length/2.0;
Rcoordinates[7][0] = 0.0; Rcoordinates[7][1] = length*sqrt(3);
Rcoordinates[8][0] = 0.0; Rcoordinates[8][1] = -length*sqrt(3);
Rcoordinates[9][0] = 2.0*length; Rcoordinates[9][1] = 0.0;
Rcoordinates[10][0] = -2.0*length; Rcoordinates[10][1] = 0.0;
Rcoordinates[11][0] = -3.0*length/2.0; Rcoordinates[11][1] = sqrt(3.0)*length/2.0;
Rcoordinates[12][0] = 3.0*length/2.0; Rcoordinates[12][1] = -sqrt(3.0)*length/2.0;
Rcoordinates[13][0] = 3.0*length/2.0; Rcoordinates[13][1] = sqrt(3.0)*length/2.0;
Rcoordinates[14][0] = -3.0*length/2.0; Rcoordinates[14][1] = -sqrt(3.0)*length/2.0;
Rcoordinates[15][0] = length; Rcoordinates[15][1] = sqrt(3.0)*length;
Rcoordinates[16][0] = -length; Rcoordinates[16][1] = -sqrt(3.0)*length;
Rcoordinates[17][0] = -length; Rcoordinates[17][1] = sqrt(3.0)*length;
Rcoordinates[18][0] = length; Rcoordinates[18][1] = -sqrt(3.0)*length;
Rcoordinates[19][0] = length + length; Rcoordinates[19][1] = sqrt(3.0)*length;
Rcoordinates[20][0] = length + 3.0*length/2.0; Rcoordinates[20][1] = sqrt(3.0)*length/2.0;
Rcoordinates[21][0] = length + 2.0*length; Rcoordinates[21][1] = 0.0;
Rcoordinates[22][0] = length + 3.0*length/2.0; Rcoordinates[22][1] = -sqrt(3.0)*length/2.0;
Rcoordinates[23][0] = length + length; Rcoordinates[23][1] = -sqrt(3.0)*length;
Rcoordinates[24][0] = 2*length + length; Rcoordinates[24][1] = sqrt(3.0)*length;
Rcoordinates[25][0] = 2*length + 3.0*length/2.0; Rcoordinates[25][1] = sqrt(3.0)*length/2.0;
Rcoordinates[26][0] = 2*length + 2.0*length; Rcoordinates[26][1] = 0.0;
Rcoordinates[27][0] = 2*length + 3.0*length/2.0; Rcoordinates[27][1] = -sqrt(3.0)*length/2.0;
Rcoordinates[28][0] = 2*length + length; Rcoordinates[28][1] = -sqrt(3.0)*length;
double shift_horz = length/2.0;
double shift_vert = sqrt(3.0)*length/2.0;
Rcoordinates[29][0] = shift_horz -length; Rcoordinates[29][1] = shift_vert  + sqrt(3.0)*length;
Rcoordinates[30][0] = shift_horz + 0.0; Rcoordinates[30][1] = shift_vert  + length*sqrt(3);
Rcoordinates[31][0] = shift_horz + length; Rcoordinates[31][1] = shift_vert  + sqrt(3.0)*length;
Rcoordinates[32][0] = shift_horz + length + length; Rcoordinates[32][1] = shift_vert + sqrt(3.0)*length;
Rcoordinates[33][0] = shift_horz + 2*length + length; Rcoordinates[33][1] = shift_vert  + sqrt(3.0)*length;
shift_horz = length/2.0;
shift_vert = -sqrt(3.0)*length/2.0;
Rcoordinates[34][0] = shift_horz -length; Rcoordinates[34][1] = shift_vert  - sqrt(3.0)*length;
Rcoordinates[35][0] = shift_horz + 0.0; Rcoordinates[35][1] = shift_vert  - length*sqrt(3);
Rcoordinates[36][0] = shift_horz + length; Rcoordinates[36][1] = shift_vert  - sqrt(3.0)*length;
Rcoordinates[37][0] = shift_horz + length + length; Rcoordinates[37][1] = shift_vert - sqrt(3.0)*length;
Rcoordinates[38][0] = shift_horz + 2*length + length; Rcoordinates[38][1] = shift_vert  - sqrt(3.0)*length;

for(int z=0; z<20; z++){
coordinates[z][0] = Rcoordinates[z][0]; coordinates[z][1] = Rcoordinates[z][1];
}

double x_shift = length/2.0;
double y_shift = length/2.0;
for(int z=20; z<40; z++){
coordinates[z][0] = Rcoordinates[z-20][0] + x_shift; coordinates[z][1] = Rcoordinates[z-20][1] + y_shift;
coordinates[z][2] = length;
}

//*************************************************************************
// Rotate structure in 3D:
//*************************************************************************

double Norm = 0.0;
double RotationNorm[theta_steps][phi_steps][omega_steps];

double XtoXTrans[theta_steps][phi_steps][omega_steps];
double XtoYTrans[theta_steps][phi_steps][omega_steps];
double XtoZTrans[theta_steps][phi_steps][omega_steps];
double YtoXTrans[theta_steps][phi_steps][omega_steps];
double YtoYTrans[theta_steps][phi_steps][omega_steps];
double YtoZTrans[theta_steps][phi_steps][omega_steps];
double ZtoXTrans[theta_steps][phi_steps][omega_steps];
double ZtoYTrans[theta_steps][phi_steps][omega_steps];
double ZtoZTrans[theta_steps][phi_steps][omega_steps];
double CenterlineProjX[theta_steps][phi_steps];
double CenterlineProjY[theta_steps][phi_steps];
double CenterlineProjZ[theta_steps][phi_steps];
double quad_place = 0.0;
double variable1, variable2;

double Proj_X, Proj_Y, Proj_Z;
double Proj_MX, Proj_MY, Proj_MZ;

double anglewt;
for(int a=0; a<theta_steps; a++){
double theta = (a*theta_step_size+theta_offset);
	for(int b=0; b<phi_steps; b++){
        double phi = (b*phi_step_size+phi_offset);
		for(int d=0; d<omega_steps; d++){//first d loop, omega loop
		double omega = (d*omega_step_size + omega_offset);

		double Trans_Matrix[3][3] = {0};
		double CenterlineProjection[3] = {0};
		RotationNorm[a][b][d] = Rotate3D(structure_tilt_angle, theta, phi, omega, Trans_Matrix, CenterlineProjection);
		XtoXTrans[a][b][d] = Trans_Matrix[0][0];
		XtoYTrans[a][b][d] = Trans_Matrix[0][1];
		XtoZTrans[a][b][d] = Trans_Matrix[0][2];
		YtoXTrans[a][b][d] = Trans_Matrix[1][0];
		YtoYTrans[a][b][d] = Trans_Matrix[1][1];
		YtoZTrans[a][b][d] = Trans_Matrix[1][2];
		ZtoXTrans[a][b][d] = Trans_Matrix[2][0];
		ZtoYTrans[a][b][d] = Trans_Matrix[2][1];
		ZtoZTrans[a][b][d] = Trans_Matrix[2][2];
		CenterlineProjX[a][b] = CenterlineProjection[0];
		CenterlineProjY[a][b] = CenterlineProjection[1];
		CenterlineProjZ[a][b] = CenterlineProjection[2];

		Norm += RotationNorm[a][b][d];
	    
            	anglewt = RotationNorm[a][b][d];

		//*************************************************************************
		// Calculate Various Projections Along Hexagonal Flake:
		//*************************************************************************

		double angle;
		double m_PointEdge[6][3] = {0.0};
		int edge = 0;
		angle = (0.0 + structure_tilt_angle)*PiNum/180.0;
		m_PointEdge[edge][0] = cos(angle)*XtoXTrans[a][b][d] + sin(angle)*YtoXTrans[a][b][d];
		m_PointEdge[edge][1] = cos(angle)*XtoYTrans[a][b][d] + sin(angle)*YtoYTrans[a][b][d];
		m_PointEdge[edge][2] = cos(angle)*XtoZTrans[a][b][d] + sin(angle)*YtoZTrans[a][b][d];
		edge = 1;
		angle = (60.0 + structure_tilt_angle)*PiNum/180.0;
		m_PointEdge[edge][0] = cos(angle)*XtoXTrans[a][b][d] + sin(angle)*YtoXTrans[a][b][d];
		m_PointEdge[edge][1] = cos(angle)*XtoYTrans[a][b][d] + sin(angle)*YtoYTrans[a][b][d];
		m_PointEdge[edge][2] = cos(angle)*XtoZTrans[a][b][d] + sin(angle)*YtoZTrans[a][b][d];
		edge = 2;
		angle = (-60.0 + structure_tilt_angle)*PiNum/180.0;
		m_PointEdge[edge][0] = cos(angle)*XtoXTrans[a][b][d] + sin(angle)*YtoXTrans[a][b][d];
		m_PointEdge[edge][1] = cos(angle)*XtoYTrans[a][b][d] + sin(angle)*YtoYTrans[a][b][d];
		m_PointEdge[edge][2] = cos(angle)*XtoZTrans[a][b][d] + sin(angle)*YtoZTrans[a][b][d];
		edge = 3;
		angle = (120.0 + structure_tilt_angle)*PiNum/180.0;
		m_PointEdge[edge][0] = cos(angle)*XtoXTrans[a][b][d] + sin(angle)*YtoXTrans[a][b][d];
		m_PointEdge[edge][1] = cos(angle)*XtoYTrans[a][b][d] + sin(angle)*YtoYTrans[a][b][d];
		m_PointEdge[edge][2] = cos(angle)*XtoZTrans[a][b][d] + sin(angle)*YtoZTrans[a][b][d];
		edge = 4;
		angle = (-120.0 + structure_tilt_angle)*PiNum/180.0;
		m_PointEdge[edge][0] = cos(angle)*XtoXTrans[a][b][d] + sin(angle)*YtoXTrans[a][b][d];
		m_PointEdge[edge][1] = cos(angle)*XtoYTrans[a][b][d] + sin(angle)*YtoYTrans[a][b][d];
		m_PointEdge[edge][2] = cos(angle)*XtoZTrans[a][b][d] + sin(angle)*YtoZTrans[a][b][d];
		edge = 5;
		angle = (180.0 + structure_tilt_angle)*PiNum/180.0;
		m_PointEdge[edge][0] = cos(angle)*XtoXTrans[a][b][d] + sin(angle)*YtoXTrans[a][b][d];
		m_PointEdge[edge][1] = cos(angle)*XtoYTrans[a][b][d] + sin(angle)*YtoYTrans[a][b][d];
		m_PointEdge[edge][2] = cos(angle)*XtoZTrans[a][b][d] + sin(angle)*YtoZTrans[a][b][d];

		//*************************************************************************
		// Pick Projection Along Hexagonal Flake clostest to X || H:
		//*************************************************************************

		int point_num = 0;
		double MagXProj = 0.0;
		double max = m_PointEdge[0][0];
		for(int pn=1; pn<6; pn++){
			if(m_PointEdge[pn][0] > max){
			max = m_PointEdge[pn][0]; point_num = pn;
			MagXProj = max;
			}//end if loop
		}//end pn loop

//*************************************************************************
//Loop over q-points:
//*************************************************************************

for(int qq=0; qq<QPoints; qq++){
double q = Data_Q[qq];//Only for C++ to SasView conversion

//*************************************************************************
// Define Cross-Sections and q-smearing:
//*************************************************************************

int qsmear_points = 5;
if(YesNo_q_smearing < 1){qsmear_points = 1;}
for(int qs=0; qs<qsmear_points; qs++){//loop for q-smearing
double QSmear_Results[2];
Calc_gaussian_smearing(Data_DeltaQ[qq], 3.0*Data_DeltaQ[qq], Data_Q[qq], Data_MeanQ[qq], qsmear_points, qs, QSmear_Results);
double Q_gaussian_wt = 1.0;
if(qsmear_points > 1){q = QSmear_Results[0]; Q_gaussian_wt = QSmear_Results[1];}

double CrossSectionVA0[4] = {0.0};
double CrossSectionVA90[4] = {0.0};
double CrossSectionPinnedSingle[4] = {0.0};

//*************************************************************************
//Define form factor scattering amplitudes:
//*************************************************************************

    double core_difference = magnetic_radius - core_radius;
    double mag_radius = core_radius + core_difference;
    double Vol = (4.0*PiNum/3.0)*pow(magnetic_radius,3);
    if(Vol == 0){Vol = 1E-10;}

    double AmpR1 = (sin(q*core_radius)/q-core_radius*cos(q*core_radius))/(q*q);
    double AmpR2 = (sin(q*(core_radius + shell1_thickness))/q-(core_radius + shell1_thickness)*cos(q*(core_radius + shell1_thickness)))/(q*q);
    double AmpR3 = (sin(q*(core_radius + shell1_thickness + shell2_thickness))/q-(core_radius + shell1_thickness + shell2_thickness)*cos(q*(core_radius + shell1_thickness + shell2_thickness)))/(q*q);
    double Amp = 4.0*PiNum*((core_sld-solvent_sld)*AmpR1 + (shell1_sld-solvent_sld)*(AmpR2-AmpR1) + (shell2_sld-solvent_sld)*(AmpR3-AmpR2));
    AmpR1 = (sin(q*mag_radius)/q-mag_radius*cos(q*mag_radius))/(q*q);
    double MAmp = 4.0*PiNum*(magcore_sld)*AmpR1;


		//*************************************************************************
		// Calculate amplitudes and cross-sections for magnetically oriented flakes:
		//*************************************************************************

		double scatt_number = 1.0*scatt_centers;
		double Norm_Volume = Vol*scatt_number;

    	    	double amp_nucl_real,amp_nucl_img,amp_mag_a_real,amp_mag_a_img,amp_mag_b_real,amp_mag_b_img,amp_mag_c_real,amp_mag_c_img;

		for(int k=0; k<2; k++){
			if(k == 0){viewing_angle = 0.0;}
			if(k == 1){viewing_angle = 90.0;}
    		double Q_X = q*cos(viewing_angle*PiNum/180.0);
    		double Q_Y = q*sin(viewing_angle*PiNum/180.0);

    	    	amp_nucl_real = 0.0;
    	    	amp_nucl_img = 0.0;
    	    	amp_mag_a_real = 0.0;
    	    	amp_mag_a_img = 0.0;
    	    	amp_mag_b_real = 0.0;
    	    	amp_mag_b_img = 0.0;
    	    	amp_mag_c_real = 0.0;
    	    	amp_mag_c_img = 0.0;

          	for(int l=0; l<scatt_centers; l++){
		Proj_X = coordinates[l][0]*XtoXTrans[a][b][d] + coordinates[l][1]*YtoXTrans[a][b][d] + coordinates[l][2]*ZtoXTrans[a][b][d];
		Proj_Y = coordinates[l][0]*XtoYTrans[a][b][d] + coordinates[l][1]*YtoYTrans[a][b][d] + coordinates[l][2]*ZtoYTrans[a][b][d];
		Proj_Z = coordinates[l][0]*XtoZTrans[a][b][d] + coordinates[l][1]*YtoZTrans[a][b][d] + coordinates[l][2]*ZtoZTrans[a][b][d];

		amp_nucl_real += Amp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_nucl_img += Amp*sin(Q_X*Proj_X + Q_Y*Proj_Y);
		if(correlation > 0){
		amp_mag_a_real += m_PointEdge[point_num][0]*MAmp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_a_img  += m_PointEdge[point_num][0]*MAmp*sin(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_b_real += m_PointEdge[point_num][1]*MAmp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_b_img  += m_PointEdge[point_num][1]*MAmp*sin(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_c_real += m_PointEdge[point_num][2]*MAmp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_c_img  += m_PointEdge[point_num][2]*MAmp*sin(Q_X*Proj_X + Q_Y*Proj_Y);

		if(MagXProj < min_parallel_projection){
		amp_mag_a_real = (MAmp/Amp)*amp_nucl_real;
		amp_mag_a_img  = (MAmp/Amp)*amp_nucl_img;
		amp_mag_b_real = 0.0;
		amp_mag_b_img  = 0.0;
		amp_mag_c_real = 0.0;
		amp_mag_c_img  = 0.0;
		}
		}//end correlation < 1 loop
		else{
		amp_mag_a_real += CenterlineProjX[a][b]*MAmp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_a_img += CenterlineProjX[a][b]*MAmp*sin(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_b_real += CenterlineProjY[a][b]*MAmp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_b_img += CenterlineProjY[a][b]*MAmp*sin(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_c_real += CenterlineProjZ[a][b]*MAmp*cos(Q_X*Proj_X + Q_Y*Proj_Y);
		amp_mag_c_img += CenterlineProjZ[a][b]*MAmp*sin(Q_X*Proj_X + Q_Y*Proj_Y);
		}//end else loop
		}//end l loop

		if(k == 0){
		viewing_angle = 0.0;
		AddToCrossSections(viewing_angle, amp_nucl_real, amp_nucl_img, amp_mag_a_real, amp_mag_a_img, amp_mag_b_real, amp_mag_b_img, amp_mag_c_real, amp_mag_c_img, anglewt, Q_gaussian_wt, Norm_Volume, CrossSectionVA0);}
		if(k == 1){
		viewing_angle = 90.0;
		AddToCrossSections(viewing_angle, amp_nucl_real, amp_nucl_img, amp_mag_a_real, amp_mag_a_img, amp_mag_b_real, amp_mag_b_img, amp_mag_c_real, amp_mag_c_img, anglewt, Q_gaussian_wt, Norm_Volume, CrossSectionVA90);}
		}//end k loop

		amp_nucl_real = Amp;
		amp_nucl_img = 0.0;
		amp_mag_a_real = CenterlineProjX[a][b]*MAmp;
		amp_mag_a_img = 0.0;
		amp_mag_b_real = CenterlineProjY[a][b]*MAmp;
		amp_mag_b_img = 0.0;
		amp_mag_c_real = CenterlineProjZ[a][b]*MAmp;
		amp_mag_c_img = 0.0;
		AddToCrossSections(0.0, amp_nucl_real, amp_nucl_img, amp_mag_a_real, amp_mag_a_img, amp_mag_b_real, amp_mag_b_img, amp_mag_c_real, amp_mag_c_img, anglewt, Q_gaussian_wt, Vol, CrossSectionPinnedSingle);

Model1[qq] += (CrossSectionVA0[0] + CrossSectionVA0[2])*volume_fraction;
Model2[qq] += (CrossSectionVA90[0] + CrossSectionVA90[2])*volume_fraction;

Model3[qq] += (pinned_frac*(CrossSectionPinnedSingle[1] + CrossSectionPinnedSingle[3]) + (1-pinned_frac)*(CrossSectionVA0[1] + CrossSectionVA0[3]))*volume_fraction;
Model4[qq] += (pinned_frac*(CrossSectionPinnedSingle[1] + CrossSectionPinnedSingle[3]) + (1-pinned_frac)*(CrossSectionVA90[1] + CrossSectionVA90[3]))*volume_fraction;

		}//end q-smearing loop
		}//end of q loop

		}//second d loop
	}//second b loop
}//second a loop

std::vector<double> m1;
std::vector<double> m2;
std::vector<double> m3;
std::vector<double> m4;

for(int qq=0; qq<QPoints; qq++){
Model1[qq] = Model1[qq]/Norm + background;
m1.push_back(Model1[qq]);
Model2[qq] = Model2[qq]/Norm + background;
m2.push_back(Model2[qq]);
Model3[qq] = Model3[qq]/Norm + SFbackground;
m3.push_back(Model3[qq]);
Model4[qq] = Model4[qq]/Norm + SFbackground;
m4.push_back(Model4[qq]);
}//end of q loop

//3. Quick plot the simulation:
//std::string title = "Hexagonally packed spheres"; char *PlotTitle = &title[0];
//std::string labelX = "q (angstroms)"; char *XLabel = &labelX[0];
//std::string labelY = "I(q)"; char *YLabel = &labelY[0];
//std::string Data1 = "NSF Horz"; char *Name1 = &Data1[0];
//std::string Data2 = "NSF Vert"; char *Name2 = &Data2[0];
//std::string Data3 = "Model NSF Horz"; char *Name3 = &Data3[0];
//std::string Data4 = "Model NSF Vert"; char *Name4 = &Data4[0];
//QuickPlotFourDataSets(4, 1, 1, QPoints, Data_Q, Data_I,  Data_I2, Model1,  Model2, PlotTitle, XLabel, YLabel, Name1, Name2, Name3, Name4);

//3. Quick plot the simulation:
//std::string title = "Hexagonally packed spheres"; char *PlotTitle = &title[0];
//std::string labelX = "q (angstroms)"; char *XLabel = &labelX[0];
//std::string labelY = "I(q)"; char *YLabel = &labelY[0];
//std::string Data1B = "SF Horz"; char *Name1B = &Data1B[0];
//std::string Data2B = "SF Vert"; char *Name2B = &Data2B[0];
//std::string Data3B = "Model SF Horz"; char *Name3B = &Data3B[0];
//std::string Data4B = "Model SF Vert"; char *Name4B = &Data4B[0];
//QuickPlotFourDataSetsB(4, 1, 1, QPoints, Data_Q, Data_I3,  Data_I4, Model3,  Model4, PlotTitle, XLabel, YLabel, Name1B, Name2B, Name3B, Name4B);
//Wait_for_keystroke();

//Write simulated data to file:
//std::string output = "Results.txt";
//char *FileToWriteOut = &output[0];
//double Data_out[5][1000] = {0};
//for(int qq=0; qq<QPoints; qq++){
//Data_out[0][qq] = Data_Q[qq];
//Data_out[1][qq] = Data_I3[qq];
//Data_out[2][qq] = Data_I4[qq];
//Data_out[3][qq] = Data_ModelSF0[qq];
//Data_out[4][qq] = Data_ModelSF90[qq];
//Write_data_noheaders(5, FileToWriteOut, Data_out, QPoints);}

json joutput = {{"m1", m1}, {"m2", m2}, {"m3", m3}, {"m4", m4}};
return joutput.dump();
}


EMSCRIPTEN_BINDINGS(my_module) {
	register_vector<double>("VectorDouble");
    emscripten::function("calculate", &calculate);
    emscripten::function("get_default_params", &get_default_params);
};





