#include "custom_functions.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>

// Function to solve the household's problem for a given starting state.
void Solve_HH_Problem(double*** Value_Function, double*** Policy_Function, int iteration, int kt0_index, int zt_index, Parameters params)
{
	// Unpack utility parameters and grids.
	double alpha = params.alpha;
	double beta = params.beta;
	double delta = params.delta;
	double* k_values = params.k_values;
	double* z_values = params.z_values;
	double** z_probs = params.z_probs;
	int number_of_k_values = params.number_of_k_values;
	int number_of_z_values = params.number_of_z_values;

	// Get starting capital and productivity values from index.
	double kt0 = k_values[kt0_index];
	double zt = z_values[zt_index];

	// Variables to store candidate optimal values.
	double v_max = -DBL_MAX;
	double kt1_optimal = 0.0;
	double new_v_max;

	for (int kt1_index = 0; kt1_index < params.number_of_k_values; kt1_index++)
	{
		// Get capital value from index.
		double kt1 = params.k_values[kt1_index];

		// Calculate value function for a given next period capital choice.
		double current_period_payoff, future_period_payoff;

		current_period_payoff = log((zt * pow(kt0, params.alpha)) + ((1 - params.delta) * kt0) - kt1);

		future_period_payoff = 0.0;
		for (int z_prob_index = 0; z_prob_index < params.number_of_z_values; z_prob_index++)
		{
			future_period_payoff += Value_Function[iteration - 1][kt1_index][z_prob_index] * params.z_probs[zt_index][z_prob_index];
		}

		new_v_max = current_period_payoff + params.beta * future_period_payoff;

		// Check if this capital choice gives the highest Value Function value.
		if (new_v_max > v_max)
		{
			// Update candidate values.
			v_max = new_v_max;
			kt1_optimal = kt1;
		}
	}

	// Update the Value Function and Policy Function with optimal values.
	Value_Function[iteration][kt0_index][zt_index] = v_max;
	Policy_Function[iteration][kt0_index][zt_index] = kt1_optimal;
}

// Function to dynamically allocate memory for a 3-dimensional array.
double*** InitializeArray3D(int size_dim_1, int size_dim_2, int size_dim_3)
{
	double*** Array3D;								// create pointer variable for start of array

	Array3D = new double** [size_dim_1];			// allocate vector of pointers for dimension 1 that will point to dimension 2 vectors

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
	Array2D = new double* [n_rows];					// create rows of pointers
	for (int i = 0; i < n_rows; i++)
	{
		Array2D[i] = new double[n_cols];			// create column of pointers for each row
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
double** ReadArrayFromCSV(std::string path)
{
	// Vector of vectors to store csv conents.
	std::vector<std::vector<double>> content_array;
	std::vector<double> row;
	std::string line, values;

	// Load file.
	std::fstream file(path, std::ios::in);
	assert(file.is_open());

	// Get data from a specific row.
	while (getline(file, line))
	{
		row.clear();

		std::stringstream row_string(line);

		// Get the individual values from each column
		while (getline(row_string, values, ','))
		{
			row.push_back(std::stod(values));		// get values in each row, and convert from string to double with stod()
		}
		content_array.push_back(row);				// add row to content array.
	}

	// Convert from std::vector to base dynamic array.
	int n_rows = content_array.size();
	int n_cols = content_array[0].size();
	double** Base_array = InitializeArray2D(n_rows, n_cols);

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			Base_array[i][j] = content_array[i][j];
		}
	}
	return Base_array;
}

// Function to read a vector from csv file (stored column-wise).
double* ReadVectorFromCSV(std::string path)
{
	// Vector of vectors to store csv conents.
	std::vector<double> content_vector;
	std::string value;

	// Load file.
	std::fstream file(path, std::ios::in);
	assert(file.is_open());

	// Go through each row
	while (getline(file, value))
	{
		// Get element, convert to double, and add to the output vector.
		content_vector.push_back(std::stod(value));
	}

	// Convert from std::vector to base dynamic array.
	int n_rows = content_vector.size();
	double* base_vector = new double[n_rows];

	for (int i = 0; i < n_rows; i++)
	{
		base_vector[i] = content_vector[i];
	}
	
	return base_vector;
}

// Alternative: read an array from file and write to a given static array.
template<size_t n_rows, size_t n_cols>
void ReadArrayFromCSV(std::string path, double(&Static_Array)[n_rows][n_cols])
{
	// Read array from file.
	double** Dynamic_Array = ReadArrayFromCSV(path);

	// Write to given array.
	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			Static_Array[i][j] = Dynamic_Array[i][j];
		}
	}
}

// Alternative: read a vector from file and write to a given static vector.
template<size_t n_rows>
void ReadVectorFromCSV(std::string path, double(&static_vector)[n_rows])
{
	// Read vector from file.
	double* dynamic_vector = ReadVectorFromCSV(path);

	// Write to given vector.
	for (int i = 0; i < n_rows; i++)
	{
		static_vector[i] = dynamic_vector[i];
	}
}