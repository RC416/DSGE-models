#include "custom_functions.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>

// Function to solve the household's problem for a given starting state.
void Solve_HH_Problem(Array3D<double>& Value_Function, Array3D<double>& Policy_Function, int iteration, int kt0_index, int zt_index, Parameters& params)
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

		// Calculate the Value Function for a given next period capital choice.
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

	// Update Value Function and Policy Function with optimal values.
	Value_Function[iteration][kt0_index][zt_index] = v_max;
	Policy_Function[iteration][kt0_index][zt_index] = kt1_optimal;
}

// Function to write a 2-dimensional array to csv file.
void WriteArrayToCSV(Array2D<double> Array, int n_rows, int n_cols, const char* file_name)
{
	std::ofstream write_output(file_name);
	assert(write_output.is_open());

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			write_output << Array[i][j] << ",";
		}
		write_output << "\n";
	}
	write_output.close();
}

// Function to read a 2-dimensional array from csv file.
Array2D<double> ReadArrayFromCSV(std::string path)
{
	// Vector of vectors to store csv conents.
	Array2D<double> content_array;
	vector<double> row;
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
	return content_array;
}

// Function to read a column vector (delimited by new lines) from csv file.
vector<double> ReadVectorFromCSV(std::string path)
{
	// Vector of vectors to store csv conents.
	vector<double> content_vector;
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
	return content_vector;
}


/*
// Function to read a row vector (delimited by commas) from csv file.
vector<double> ReadVectorFromCSV(std::string path)
{
	// Vector of vectors to store csv conents.
	vector<double> content_vector;
	std::string value;

	// Load file.
	std::fstream file(path, std::ios::in);
	assert(file.is_open());

	// Go through the row of data, stopping at each ',' to get the value.
	while (getline(file, value, ','))
	{
		// Get element, convert to double, and add to the output vector.
		content_vector.push_back(std::stod(value));
	}
	return content_vector;
}
*/

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