#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "cfi.h"
#include "Matrix.h"
#include "Coord.h"
#include "Table.C"
#include "geo_etc.h"

using namespace std;

/*
 * beamformer input_file rawdata_file output_file
 *
 * input_file:
 *   Source x, y, z
 *   Receiver x, y, z
 *   Receiver length
 *   Receiver heading  // in degree from utm_x in ccw direction
 *                     // converted to radian inside the code
 *   range_resolution
 *   beamform_freq
 *   two_way_range_of_interest
 *   angle_of_interest  // in degree from receiver, measured from utm_x in 
 *                      // ccw direction
 *                      // converted to radian inside the code
 */

int main(int argc, char *argv[])
{

	if (argc != 4)
	{
		cerr << "error in argument" << endl;
 		cerr << "beamfomer input_file rawdata_file output_file" << endl;
		exit(1);
	}

	CartesianCoord SourcePosition;
	_Array Receiver;
	float range_resolution;
	float beamform_freq;
	float range_start, range_end;
	float angle_start, angle_end;

	string rawdata_filename, output_filename;

	ofstream log_file("beamformer.log", ios::out);

	////////////////////////////////////////////////////////////////////////
	// Input
	////////////////////////////////////////////////////////////////////////
	FILE * fin = fopen(argv[1], "rt");
	char Buffer[1000];
	float temp_float, temp_x, temp_y, temp_z;
	int temp_int;

	// source
	CommentFileIn(fin, temp_x);
	CommentFileIn(fin, temp_y);
	CommentFileIn(fin, temp_z);
	SourcePosition.SetCoord(temp_x, temp_y, temp_z);
	log_file << "SourcePosition: " << SourcePosition.GetX() << ", " << 
			SourcePosition.GetY() << ", " << SourcePosition.GetZ() << endl;

	// receiver
	CommentFileIn(fin, temp_x);
	CommentFileIn(fin, temp_y);
	CommentFileIn(fin, temp_z);
	Receiver.Position.SetCoord(temp_x, temp_y, temp_z);
	CommentFileIn(fin, temp_float);
	Receiver.Length = temp_float;
	CommentFileIn(fin, temp_float);
	Receiver.Heading = temp_float * M_PI / 180.;
	log_file << "ReceiverPosition: " << Receiver.Position << endl;
	log_file << "ReceiverLength: " << Receiver.Length << endl;
	log_file << "ReceiverHeading: " << Receiver.Heading << endl;

	// range_resolution
	CommentFileIn(fin, range_resolution);
	log_file << "range_resolution: " << range_resolution << endl;

	// beamform frequency
	CommentFileIn(fin, beamform_freq);
	log_file << "Beamform freq: " << beamform_freq << endl;

	// range
	CommentFileIn(fin, range_start);
	CommentFileIn(fin, range_end);
	log_file << "Range start, end: " 
		<< range_start << ", " << range_end << endl;

	// angle
	CommentFileIn(fin, angle_start);
	CommentFileIn(fin, angle_end);
	angle_start *= (M_PI / 180.);
	angle_end *= (M_PI / 180.);
	log_file << "Angle start, end (deg.): " << 
			angle_start * 180. / M_PI << ", " << 
			angle_end  * 180. / M_PI << endl;

	// rawdata_file
	rawdata_filename = string(argv[2]);

	// output_file
	output_filename = string(argv[3]);

	//////////////////////////////////////////////////////////////
	// Read RawData;
	//////////////////////////////////////////////////////////////

	_ARR RawDataArr;
	cout << "Reading data:" << endl; cout.flush();
	RawDataArr.ReadData(rawdata_filename);
	cout << "done" << endl; cout.flush();

	//////////////////////////////////////////////////////////////
	// Set Elliptic Map Dimension
	//////////////////////////////////////////////////////////////
	
	int r_size, theta_size;
	int grid_r_min, grid_r_max;
	int grid_mili_theta_min, grid_mili_theta_max;
	int grid_r_inc, grid_mili_theta_inc;

	// ## points per range resolution
	grid_r_inc = int(range_resolution / 2. /5.);

	float lambda;
	const float sound_speed = 1500.;
	float K_Lambda_over_L;
	float Endfire_Resolution;
	lambda = sound_speed / beamform_freq;
	K_Lambda_over_L = 1.3 * lambda / Receiver.Length;
	Endfire_Resolution = 2.16 * 1.3 * sqrt(lambda / Receiver.Length);

	// at least ## points per cross range resolution
	grid_mili_theta_inc = 40;
	//	int(K_Lambda_over_L * 180. / M_PI * 1000 / 2. / 10.);

	grid_r_min = int(range_start - range_resolution);
	grid_r_max = int(range_end + range_resolution) + 1;

	grid_mili_theta_min = int(angle_start * 180. / M_PI * 1000 
			- Endfire_Resolution * 180. / M_PI * 1000 / 2.);
	grid_mili_theta_max = int(angle_end * 180. / M_PI * 1000 
			+ Endfire_Resolution * 180. / M_PI * 1000 / 2.) + 1;

	_EllipticARR EllipticData;

	EllipticData.grid_r_inc = grid_r_inc;
	EllipticData.grid_mili_theta_inc = grid_mili_theta_inc;
	EllipticData.grid_r_min = grid_r_min;
	EllipticData.grid_r_max = grid_r_max;
	EllipticData.grid_mili_theta_min = grid_mili_theta_min;
	EllipticData.grid_mili_theta_max = grid_mili_theta_max;

	cout << "Memory allocation..."; cout.flush();
	EllipticData.SetDimension();
	cout << "done" << endl; cout.flush();

	//////////////////////////////////////////////////////////////
	// Set Elliptic Coordinates
	//////////////////////////////////////////////////////////////

	EllipticData.SourcePosition.SetX(SourcePosition.GetX());
	EllipticData.SourcePosition.SetY(SourcePosition.GetY());
	EllipticData.ReceiverPosition.SetX(Receiver.Position.GetX());
	EllipticData.ReceiverPosition.SetY(Receiver.Position.GetY());

	EllipticData.CenterPosition = 
		(EllipticData.SourcePosition + EllipticData.ReceiverPosition) / 2.;

	//////////////////////////////////////////////////////////////
	// Import raw data into elliptic coord
	//////////////////////////////////////////////////////////////

	cout << "Importing data... "; cout.flush();
	Table<PolarCoord> PolarFromReceiver;
	EllipticData.ImportArr(RawDataArr, PolarFromReceiver);
	cout << "done" << endl; cout.flush();

///	EllipticData.WriteData("temp_elliptic.arr");

	//////////////////////////////////////////////////////////////
	// Beamform
	//////////////////////////////////////////////////////////////
	
	_EllipticARR BeamformedElliptic;

	cout << "Memory allocation..."; cout.flush();
	BeamformedElliptic.Data.Dimension(EllipticData.r_size, 
										EllipticData.theta_size);
	cout << "done" << endl; cout.flush();
	
	cout << "Beamforming... ";
	BeamformedElliptic = 
		EllipticData.Beamform(PolarFromReceiver, Receiver, range_resolution, 
				Endfire_Resolution, beamform_freq);
	cout << "done" << endl;

///	BeamformedElliptic.WriteData("temp_elliptic.arr");

	//////////////////////////////////////////////////////////////
	// Export beamfored field into Rectangular ARR
	//////////////////////////////////////////////////////////////
	_ARR BeamformedArr;
	BeamformedArr.footer_found = 1;
	BeamformedArr.grid_inc = RawDataArr.grid_inc;
	BeamformedArr.grid_xmin = RawDataArr.grid_xmin;
	BeamformedArr.grid_xmax = RawDataArr.grid_xmax;
	BeamformedArr.grid_ymin = RawDataArr.grid_ymin;
	BeamformedArr.grid_ymax = RawDataArr.grid_ymax;
	BeamformedArr.Data.Dimension(
			RawDataArr.Data.GetRow(), RawDataArr.Data.GetColumn());
	BeamformedArr.Data.Fill(0.);

	cout << "Exporting data... ";
	BeamformedElliptic.Export2Arr(BeamformedArr);
	cout << "done" << endl;

	BeamformedArr.WriteData(output_filename);


	log_file.close();
	return 0;
}
