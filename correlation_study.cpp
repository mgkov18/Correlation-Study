/*
	Programmer: Mike Kovacevich
	Email: kovacevich.9@osu.edu
	Correlation Study using data from IceCube 7 year survey

	 Revision History:
	 4/23/17 Created functions to read in data from dat files in 2-D vectors
	 4/27/17 Tried to implement a new approach using strings to make the 2-D vectors
	 4/28/17 Implemented a new way of making 2-D vectors, the way that I first tried on 4/23/17
	 		 lead to multiple bugs and the 2-D vectors not being created properly
	 4/29/17 Added the gsl routine to help calculate the surface area of the band and the surface
	 		  area of the patch in the band 
	 5/3/17  Sorting through 2-D vectors in order to later group events
	 5/12/17 Created a matrix that contains all the neutrino events
	 5/16/17 Implementing a function that calculates the distance between a neutrino and cosmic ray,
	 		 using the formula for arc length and angle is in radians
	 		 Also removed the gsl integration routines
	 5/23/17 Commented out certain lines of code, deleted the (unnecessay?) 9th coumn of data in the IC86-year dat files
			 Also deleted unnecessary columns in the dat file that contains the cosmic ray data
	 5/26/17 deleted portions of the events_window function and attempted to load neutrino data from
	 		 a matrix with all the neutrinos to a matrix that only contains the R.A. and declination
	 ***** Matrix = 2-D Vector ******
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string> 
#include <gsl/gsl_integration.h>

using namespace std;

void events_window(vector< vector<double> > &neutrinos, vector< vector<double> > &cosmic_rays,
				   double radius);

int main()
{
	  const double total_neutrino_events = 898281; //total number of neutrino events
	  const double total_cosmic_ray_events = 1233487; //total number of cosmic_ray_events
	  const double earth_radius = 6.37e6; 
	  const double cosmic_columns = 9; //the number of columns in the cosmic ray dat file
	  const double neutrino_columns = 8; //number of columns in the neutrino dat files
  	
	//Creating a 2-D vector that reads in cosmic ray data from a dat file
	vector<vector<double> > cosmic_ray_events(total_cosmic_ray_events, vector<double>(9,0));
	ifstream cosmic_ray_data("events_comp-h4a_10.00-10000.00PeV_zen37.00.dat", ios::in);   
	//checking that the dat file exists/open so that the data can be streamed as input    
	while(cosmic_ray_data.good())    
	{
		for(int i = 0; i < 1233487; i++) //1233487 is the number of rows in the dat file
		{
			for(int j = 0; j < cosmic_columns; j++) 
			{	
					cosmic_ray_data >> cosmic_ray_events[i][j]; //reading in data for 2-D vector
			}
		}
	}

    /*
	The header information for the cosmic ray data is listed on the line below:
	1)MJD  2)second within MJD  3) S125[VEM]  4)Zenith angle (rad)  5)Azimuthal angle (rad)
	6)R.A. (rad)  7)Declination (rad)
    */

	cout << "Printing out values from the cosmic rays matrix: " << endl << endl;
		//checking the cosmic ray matrix is made correctly		 
		for(int i = 0; i < 9; i++)
		{
			for(int j = 0; j < 9; j++) 
			{
				cout << cosmic_ray_events[i][j] << "    ";
			}
			cout << endl;
		}
     
	cout << endl << "Comic ray values have been printed out" << endl << endl;

   	
	cosmic_ray_data.close();	

	//Combining all the neutrino data into one file
	ofstream out;
	out.open("all_neutrino_events.dat");

	//Line 82 is commented out because it distrupts the creation of the 2-D vector containing
	//all the neutrino events but Line 82 contains the header information
	//out << "Run    Event     R.A.     Dec      log(E)    Sigma     Time      sinDec" << endl;

	//declaring a 36900x8 vector, dat file has 36900 rows
	vector<vector<double> > IC40_events (36900, vector<double>(neutrino_columns,0));
	ifstream IC40_events_data("IC40_exp.dat", ios::in);
	//the while loop makes sure that the dat file is open/exists
	while(IC40_events_data.good())
	{
		for(int i = 0; i < 36900; i++) //This dat file contained 36900 rows 
		{	
			for (int j = 0; j < neutrino_columns; j++) 
			{
				IC40_events_data >> IC40_events[i][j];  //building the 2-D vector
			}
		}
	}
	IC40_events_data.close();
		//writing the 2-D vector to a dat file with all the neutrino events
		for(int i = 0; i < 36900; i++)  //was 36900
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC40_events[i][j] << "      ";
			}
			out << endl;
		}

	//dat file contained 107012 rows, declared a 107012x8 vector
	vector<vector<double> > IC59_events(107011, vector<double>(neutrino_columns,0));
	ifstream IC59_events_data("IC59_exp.dat", ios::in);
	//checking the dat file exists and input can be streamed from it
	while(IC59_events_data.good())
	{
		for(int i = 0; i < 107011; i++)
		{
			for (int j = 0; j < neutrino_columns; j++)
			{
				IC59_events_data >> IC59_events[i][j];
			}
		}
	}
	IC59_events_data.close();

		for(int i = 0; i < 107011; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC59_events[i][j] << "      ";
			}
			out << endl;
		}


	//delcaring a 93134x8 vector, dat file contained 93134 rows
	vector<vector<double> > IC79b_events(93133, vector<double>(neutrino_columns,0));
	ifstream IC79b_events_data("IC79b_exp.dat", ios::in);
	//checking that the data in the dat file can be streamed as input/the file exists
	while(IC79b_events_data.good())
	{
		for(int i = 0; i < 93133; i++)
		{
			for (int j = 0; j < neutrino_columns; j++)
			{
				IC79b_events_data >> IC79b_events[i][j]; //building the vector
			}
		}
	}
	IC79b_events_data.close();

		for(int i = 0; i < 93133; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC79b_events[i][j] << "      ";
			}
			out << endl;
		}

	//declaring a 109741x8 vector, dat file contained 109741 rows
	vector<vector<double> > IC79_events(109740, vector<double>(neutrino_columns,0));
	ifstream IC79_events_data("IC79_exp.dat", ios::in);
	//checking existence of the file
	while(IC79_events_data.good())
	{
		for(int i = 0; i < 109740; i++)
		{
			for (int j = 0; j < neutrino_columns; j++)
			{
				IC79_events_data >> IC79_events[i][j]; //building the vector
			}
		}
	}
	IC79_events_data.close();

		for(int i = 0; i < 109740; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC79_events[i][j] << "      ";
			}
			out << endl;
		}

	//declaring a 105301x8 vector, dat file contained 105301 rows
	vector<vector<double> > IC86_2012_events(105300, vector<double>(9,0));
	ifstream IC86_2012_events_data("IC86-2012_exp.dat", ios::in);
	//checking existence of the file
	while(IC86_2012_events_data.good())
	{
		for(int i = 0; i < 105300; i++) 
		{
			for (int j = 0; j < 8; j++)
			{
				IC86_2012_events_data >> IC86_2012_events[i][j]; //Adding in the elements of the vector
			}
		}
	}
	IC86_2012_events_data.close();

		for(int i = 0; i < 105300; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC86_2012_events[i][j] << "      ";
			}
			out << endl;
		}		

	//declaring a 114835x8 vector, dat file contained 114835 rows
	//this dat file contained an extra column of data, although it does not impact the calculations
	//that I will do
	vector<vector<double> > IC86_2013_events(114834,vector<double>(9,0));
	ifstream IC86_2013_events_data("IC86-2013_exp.dat", ios::in);
	//checking existence of the file
	while(IC86_2013_events_data.good())
	{
		for(int i = 0; i < 114834; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				IC86_2013_events_data >> IC86_2013_events[i][j]; //building the 2-D vector
			}
		}
	}
	IC86_2013_events_data.close();

		for(int i = 0; i < 114834; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC86_2013_events[i][j] << "      ";
			}
			out << endl;
		}	

	//declaring a 118457x8 vector, dat file contained 118457
	vector<vector<double> > IC86_2014_events(118456,vector<double>(9,0));
	ifstream IC86_2014_events_data("IC86-2014_exp.dat", ios::in);
	//checking existence of the file
	while(IC86_2014_events_data.good())
	{
		for(int i = 0; i < 118456; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				IC86_2014_events_data >> IC86_2014_events[i][j]; //building the 2-D vector
			}
		}
	}
	IC86_2014_events_data.close();

		for(int i = 0; i < 118456; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC86_2014_events[i][j] << "      ";
			}
			out << endl;
		}

	//declaring the 76664x8 vector, dat file contained 76664 rows
	vector<vector<double> > IC86_2015_events(76663,vector<double>(9,0));
	ifstream IC86_2015_events_data("IC86-2015_exp.dat", ios::in);
	//checking the file is open/exists and can be streamed as input
	while(IC86_2015_events_data.good())
	{
		for(int i = 0; i < 76663; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				IC86_2015_events_data >> IC86_2015_events[i][j]; //building the 2-D vector
			}
		}
	}
	IC86_2015_events_data.close();

		for(int i = 0; i < 76663; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC86_2015_events[i][j] << "      ";
			}
			out << endl;
		}

	//declaring a 2-D vector of size 136245x8, dat file contained 136245 rows
	vector<vector<double> > IC86_events(136244,vector<double>(neutrino_columns,0));  
	ifstream IC86_events_data("IC86_exp.dat", ios::in);
	//checking the dat file is open and input can be streamed into a vector
	while(IC86_events_data.good())
	{
		for(int i = 0; i < 136244; i++)
		{
			for (int j = 0; j < neutrino_columns; j++)
			{
				IC86_events_data >> IC86_events[i][j]; //building the vector
			}
		}
	}
	IC86_events_data.close();

		for(int i = 0; i < 136244; i++)
		{
			for(int j = 0; j < 8; j++)
			{
				out << IC86_events[i][j] << "     ";
			}
			out << endl;
		}

  	//checks that the dat files were wriiten to with the right data
  cout << "All neutrino events written to neutrino_events.dat" << endl;
  out.close();

    vector<vector<double> > all_neutrinos(898281,vector<double>(8,0)); 
    vector<vector<double> > sorted_neutrinos
  
	ifstream all_neutrino_events("all_neutrino_events.dat", ios::in);
	//checking the dat file is open and input can be streamed into a matrix
	while(all_neutrino_events.good())
	{
		for(int i = 0; i < 898281; i++)  //was i < 346784
		{
			for (int j = 0; j < neutrino_columns; j++)
			{
				if(j == 2)
				{
					all_neutrino_events >> all_neutrinos[i][2]; //building the vector
				}
				else if(j == 3)
				{
					all_neutrino_events >> all_neutrinos[i][3];
				}
			}
		}
	}
	all_neutrino_events.close();

	cout << "R.A. and declination have been sorted" << endl;

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			cout << all_neutrinos[i][j] << "   ";
		}
		cout << endl;
	}

	cout << endl << "All neutrino events and R.A. and declination are now loaded into in a matrix (2-D vector)" << endl << endl;

	//events_window(all_neutrinos, cosmic_ray_events, earth_radius);
	
	return(0);
}

//This function calculates the distance between two points, mainly a neutrino and cosmic ray and 
//this function also checks if the cosmic ray is in the window
void events_window(vector< vector<double> > &neutrinos, vector< vector<double> > &cosmic_rays,
					 double radius)
{
	//selecting certain elements from the dat file containing the cosmic ray data
	ifstream cosmic_ray_angle_out;
	cosmic_ray_angle_out.open ("events_comp-h4a_10.00-10000.00PeV_zen37.00.dat");

	for(int i = 0; i < 1233487; i++)
	{
		for(int j = 0; j < 9; j++)
		{
			if(j = 5)  //selects the 5th element since this contains the Right Ascension
				cosmic_ray_angle_out >> cosmic_rays[i][j];
			else if(j = 6)  //selects the 6th element since this contains the Declination
				cosmic_ray_angle_out >> cosmic_rays[i][j];
		}
	}

	cout << "Checking the cosmic rays matrix and selecting only certain values" << endl;
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			cout << cosmic_rays[i][j];
		}
		cout << endl;
	}

	ofstream event_count;
	event_count.open("events_count.dat");

	//calculating the number of neutrinos in a band and window by narrowing down the number of 
	//events that occur at differing angles

	
	for(int i = 0; i < 898281; i++)
	{
		int count = 0;
		for(int z = 0; z < 3; z++) //should be less than 1233487
		{
			event_count << neutrinos[i][3] << endl;
		}

		for(int j = 0; j < 10; j++)
		{
			for(int k = 0; k < 3; k++)
			{
				if(k == 2)
				{
					double arc_length = radius * (neutrinos[i][2] - cosmic_rays[j][2]);
						if(arc_length <= 5.) //.05 is an arbitrary number, real number will be obtained
						{                     //from the Monte Carlo simulation
							event_count << count++ << "    "  << cosmic_rays[j][2] << endl; 
						}
						else{}
				}
				else{}
			}	
		}
		event_count << count << endl;
	}

	
	event_count.close();
	cosmic_ray_angle_out.close();
}