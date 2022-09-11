// Base Model implementation in C++.

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <cfloat>
#include <vector>

// Declarations.
using std::vector;
void WriteArrayToCSV(vector<vector<double>> Array2D, const char* file_name);

int main()
{
	// Assign parameter values.
	double alpha = 0.400;
	double beta = 0.987;
	double delta = 1.000;
	int const number_of_iterations = 1000;

	// Calculate the steady-state level of capital.
	double k_steady = pow(((1 - beta * (1 - delta)) / (alpha * beta)), (1 / (alpha - 1)));  // pow(base, exponent) = base^exponent

	// Create grid of capital values around the steady-state (+/- 50%).
	int const number_of_k_values = 201;
	double k_low_pct = 0.50;
	double k_high_pct = 1.50;
	double k_values[number_of_k_values];

	for (int i = 0; i < number_of_k_values; i++)
	{
		k_values[i] = k_steady * (k_low_pct + ((double(i) / (double(number_of_k_values) - 1)) * (k_high_pct - k_low_pct)));
	}

	// Initialize the Value Function and Policy Function (as arrays).
	vector<vector<double>> Value_Function(number_of_iterations, vector<double>(number_of_k_values, 0));
	vector<vector<double>> Policy_Function(number_of_iterations, vector<double>(number_of_k_values, 0));

	// Solve the household's problem for each possible starting state.
	for (int iteration = 1; iteration < number_of_iterations; iteration++)
	{
		for (int kt0_ind = 0; kt0_ind < number_of_k_values; kt0_ind++)
		{
			// Variables to store candidate optimal values.
			double v_max = -DBL_MAX;
			double kt1_optimal = 0.0;
			double new_value_function_value;

			for (int kt1_ind = 0; kt1_ind < number_of_k_values; kt1_ind++)
			{
				// Calculate the Value Function for given starting capital and next period capital choice.
				new_value_function_value = log(pow(k_values[kt0_ind], alpha) - k_values[kt1_ind] + (1 - delta) * k_values[kt0_ind]) +
					+beta * Value_Function[iteration - 1][kt1_ind];

				// Check if this capital choice gives highest Value Function value.
				if (new_value_function_value > v_max)
				{
					// Update candidate values.
					v_max = new_value_function_value;
					kt1_optimal = k_values[kt1_ind];
				}
			}

			// Update the Value Function and Policy Function with optimal values.
			Value_Function[iteration][kt0_ind] = v_max;
			Policy_Function[iteration][kt0_ind] = kt1_optimal;
		}
	}

	// Write Value Function and Policy Function to csv files.
	WriteArrayToCSV(Value_Function, "Value_Function.csv");
	WriteArrayToCSV(Policy_Function, "Policy_Function.csv");

	// Display the first and last 5 values of the Value Function and Policy Function.
	for (int i = 0; i < 5; i++)
	{
		if (i == 0) { std::cout << "Value Function V(k):" << "\n"; };
		std::cout << "V(" << k_values[i] << ") = " << Value_Function[number_of_iterations - 1][i] << "\n";
	}
	std::cout << "..." << "\n";
	for (int i = number_of_k_values - 5; i < number_of_k_values; i++)
	{
		std::cout << "V(" << k_values[i] << ") = " << Value_Function[number_of_iterations - 1][i] << "\n";
	}

	for (int i = 0; i < 5; i++)
	{
		if (i == 0) { std::cout << "\n" << "Policy Function g(k):" << "\n"; };
		std::cout << "g(" << k_values[i] << ") = " << Policy_Function[number_of_iterations - 1][i] << "\n";
	}
	std::cout << "..." << "\n";
	for (int i = number_of_k_values - 5; i < number_of_k_values; i++)
	{
		std::cout << "g(" << k_values[i] << ") = " << Policy_Function[number_of_iterations - 1][i] << "\n";
	}

	// Leave window open after program terminates.
	std::cin.get();
	return 0;
}

// Function to write a 2-dimensional array to csv file.
void WriteArrayToCSV(vector<vector<double>> Array2D, const char* file_name)
{
	int n_rows = Array2D.size();
	int n_cols = Array2D[0].size();

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
