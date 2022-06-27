// Base Model implementation in C++

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

double** InitiateArray_2D(int n_rows, int n_cols);
void DeleteArray_2D(double** Array_2D, int n_rows, int n_cols);
void WriteArrayToCSV(double** Array_2D, int n_rows, int n_cols, const char* file_name);

int main()
{
	// Assign initial values.
	const double alpha = 0.400;
	const double beta = 0.987;
	const double delta = 1.000;
	const int number_of_iterations = 1000;

	// Calculate steady-state level of capital.
	double k_steady = pow(((1 - beta * (1 - delta)) / (alpha * beta)), (1 / (alpha - 1)));

	// Create range of capital values around steady-state.
	const int number_of_k_values = 201;
	const double k_low_pct = 0.50;
	const double k_high_pct = 1.50;
	double k_values[number_of_k_values];

	for (int i = 0; i < number_of_k_values; i++)
	{
		k_values[i] = k_low_pct*k_steady + (double(i) / (number_of_k_values-1)) * ((k_high_pct - k_low_pct)*k_steady);
	}

	// Initialize Value Function and Policy Function (as arrays).
	double** Value_Function = InitiateArray_2D(number_of_iterations, number_of_k_values);
	double** Policy_Function = InitiateArray_2D(number_of_iterations, number_of_k_values);
	for (int i=0;i<number_of_k_values;i++)
	{
		Value_Function[0][i] = 0.0; // assign value of 0 to first Value Function iteration
	}

	// Perform value function iteration.
	for (int iteration = 1; iteration < number_of_iterations; iteration++)
	{	
		// double v_max, kt1_optimal, new_value_function_value;

		for (int kt0_ind = 0; kt0_ind < number_of_k_values; kt0_ind++)		// for each level of starting capital...
		{	
			// Variables to store candidate optimal values.
			double v_max = -DBL_MAX;
			double kt1_optimal = 0.0;
			double new_value_function_value;

			for (int kt1_ind = 0; kt1_ind < number_of_k_values; kt1_ind++)	// ... check all next period capital choices
			{
				// Calculate Value Function for given starting capital and next period capital choice.
				new_value_function_value = log(pow(k_values[kt0_ind], alpha) - k_values[kt1_ind] + (1 - delta) * k_values[kt0_ind] ) +
					+ beta * Value_Function[iteration-1][kt1_ind];

				// Check if this capital choice gives highest Value Function value.
				if (new_value_function_value > v_max)
				{
					// Update candidate values.
					v_max = new_value_function_value;
					kt1_optimal = k_values[kt1_ind];
				}
			}

			// Update Value Function and Policy Function with optimal values.
			Value_Function[iteration][kt0_ind] = v_max;
			Policy_Function[iteration][kt0_ind] = kt1_optimal;
		}
	}

	// Write Value Function and Policy Function to csv files.
	WriteArrayToCSV(Value_Function, number_of_iterations, number_of_k_values, "Value_Function.csv");
	WriteArrayToCSV(Policy_Function, number_of_iterations, number_of_k_values, "Policy_Function.csv");


	// Display the first and last 5 values of the Value Function and Policy Function.
	for (int i = 0; i < 5; i++)
	{
		if (i == 0) { std::cout << "Value Function V(k):" << "\n"; };
		std::cout << "V(" << k_values[i] << ") = " << Value_Function[number_of_iterations - 1][i] << "\n";
	}
	std::cout << "..." << "\n";
	for (int i = number_of_k_values-5; i < number_of_k_values; i++)
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

	// Remove arrays from memory.
	DeleteArray_2D(Value_Function, number_of_iterations, number_of_k_values);
	DeleteArray_2D(Policy_Function, number_of_iterations, number_of_k_values);
	
	std::cin.get();
	return 0;
}


// Function to dynamically allocate memory for a 2-dimensional array.
double** InitiateArray_2D(int n_rows, int n_cols)
{
	double** Array_2D;
	Array_2D = new double* [n_rows];		// create rows of pointers
	for (int i = 0; i < n_rows; i++)
	{
		Array_2D[i] = new double[n_cols];	// create column of pointers for each row
	}

	return Array_2D;
}

// Function to clear the memory allocated for a 2-dimensional array.
void DeleteArray_2D(double** Array_2D, int n_rows, int n_cols)
{
	for (int i = 0; i < n_rows; i++)
	{
		delete[] Array_2D[i];
	}
	delete[] Array_2D;
}

// Function to write a 2-dimensional array to csv file.
void WriteArrayToCSV(double** Array_2D, int n_rows, int n_cols, const char* file_name)
{
	std::ofstream write_output(file_name);
	assert(write_output.is_open());

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			write_output << Array_2D[i][j] << ",";
		}
		write_output << "\n";
	}
	write_output.close();
}
