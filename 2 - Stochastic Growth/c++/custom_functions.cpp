#include "custom_functions.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>

// Function to solve Value Function for given starting states.
void Iterate_Value_Function(double*** Value_Function, double*** Policy_Function, int iteration, int kt0_index, int zt_index, Parameters params)
{
	// Get starting capital and productivity values from index.
	double kt0 = params.k_values[kt0_index];
	double zt = params.z_values[zt_index];

	// Variables to store candidate optimal values.
	double v_max = -DBL_MAX;
	double kt1_optimal = 0.0;
	double new_v_max;

	for (int kt1_index = 0; kt1_index < params.number_of_k_values; kt1_index++)
	{
		// Get capital value from index.
		double kt1 = params.k_values[kt1_index];

		// Calculate value function for given choice of next period capital.
		double current_period_payoff, future_period_payoff;

		current_period_payoff = log((zt * pow(kt0, params.alpha)) + ((1 - params.delta) * kt0) - kt1);

		future_period_payoff = 0.0;
		for (int z_prob_index = 0; z_prob_index < params.number_of_z_values; z_prob_index++)
		{
			future_period_payoff += Value_Function[iteration - 1][kt1_index][z_prob_index] * params.z_probs[zt_index][z_prob_index];
		}

		new_v_max = current_period_payoff + params.beta * future_period_payoff;

		// Check if this capital choice gives highest Value Function value.
		if (new_v_max > v_max)
		{
			// Update candidate values.
			v_max = new_v_max;
			kt1_optimal = kt1;
		}
	}

	// Update Value Function and Policy Function with optimal values.
	Value_Function[iteration][kt0_index][zt_index] = v_max;
	Policy_Function[iteration][kt0_index][zt_index] = kt1_optimal;
}

// Function to dynamically allocate memory for a 3-dimensional array.
double*** InitializeArray3D(int size_dim_1, int size_dim_2, int size_dim_3)
{
	double*** Array3D;								// create pointer variable for start of array

	Array3D = new double** [size_dim_1];		// allocate vector of pointers for dimension 1 that will point to dimension 2 vectors

	for (int i = 0; i < size_dim_1; i++)
	{
		Array3D[i] = new double* [size_dim_2];		// allocate vector of pointers for dimension 2 that will point to dimension 3 vectors

		for (int j = 0; j < size_dim_2; j++)
		{
			Array3D[i][j] = new double[size_dim_3]; // allocate vector of doubles for dimension 3
		}
	}
	return Array3D;
}

// Function to clear the memory allocated for a 3-dimensional array.
void DeleteArray3D(double*** Array3D, int size_dim_1, int size_dim_2, int size_dim_3)
{
	for (int i = 0; i < size_dim_1; i++)
	{
		for (int j = 0; j < size_dim_2; j++)
		{
			delete[] Array3D[i][j];					// delete each dimension 3 vector of data
		}
		delete[] Array3D[i];						// delete each dimension 2 vector of pointers
	}
	delete[] Array3D;								// delete pointer to dimension 1 vector of pointers
}

// Function to dynamically allocate memory for a 2-dimensional array.
double** InitializeArray2D(int n_rows, int n_cols)
{
	double** Array2D;
	Array2D = new double* [n_rows];		// create rows of pointers
	for (int i = 0; i < n_rows; i++)
	{
		Array2D[i] = new double[n_cols];	// create column of pointers for each row
	}

	return Array2D;
}

// Function to clear the memory allocated for a 2-dimensional array.
void DeleteArray2D(double** Array2D, int n_rows, int n_cols)
{
	for (int i = 0; i < n_rows; i++)
	{
		delete[] Array2D[i];
	}
	delete[] Array2D;
}

// Function to write a 2-dimensional array to csv file.
void WriteArrayToCSV(double** Array2D, int n_rows, int n_cols, const char* file_name)
{
	std::ofstream write_output(file_name);
	assert(write_output.is_open());

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			write_output << Array2D[i][j] << ",";
		}
		write_output << "\n";
	}
	write_output.close();
}

// Function to read a 2-dimensional array from csv file.
// Uses std::vector so that csv array dimensions do not need to be specified.
std::vector<std::vector<double>> ReadArrayFromCSV(std::string path)
{
	// Vector of vectors to store csv conents.
	std::vector<std::vector<double>> Content;
	std::vector<double> row;
	std::string line, values;

	// Load file.
	std::fstream file(path, std::ios::in);
	assert(file.is_open());

	while (getline(file, line))
	{
		row.clear();

		std::stringstream row_string(line);

		while (getline(row_string, values, ','))
		{
			row.push_back(std::stod(values));		// get values in each row, and convert from string to double with stod()
		}
		Content.push_back(row);							// add row to content array.
	}

	/*
	// Print array.
	for (int i = 0; i < Content.size(); i++)
	{
		for (int j = 0; j < Content[0].size(); j++)
		{
			std::cout << Content[i][j] << " ";
		}
		std::cout << "\n";
	}
	*/

	return Content;
}