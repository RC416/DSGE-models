#include "custom_functions.h"
#include <iostream>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>

/*
 -----------------------------------------------------------------------------------------------------
 Function 1 - Solve Household Problem.
 -----------------------------------------------------------------------------------------------------
Function to find the next iteration of the Value Function and Policy Function
by solving the household's problem for a given starting state.
Outputting multiple values is handled by passing references to output variables.

Input:
	- Current Value Function
	- Market price for bond
	- Starting wealth and endowment level
	- References to output variables (V, g, g_index)

Output:
	- Next value function value
	- Optimal choice of savings/borrowing for next period
 */

void Solve_HH_Problem(double q, int a_start_index, int e_start_index, Array2D<double>& Value_Function, 
					  Parameters& params, double& V, double& g, int& g_index)
{
	// Unpack utility parameters and grids.
	double sigma = params.sigma;
	double beta = params.beta;
	vector<double> a_grid = params.a_grid;
	vector<double> e_grid = params.e_grid;
	Array2D<double> e_probs = params.e_probs;
	int number_of_a_values = params.number_of_a_values;
	int number_of_e_values = params.number_of_e_values;

	// Get value of state variables.
	double a_start = a_grid[a_start_index];
	double e_start = e_grid[e_start_index];

	// Variables to store candidate values.
	double v_max = -DBL_MAX;
	double new_v_max;
	double a_next_optimal = 0.0;
	int a_next_optimal_index = 0;
	double current_period_payoff, future_period_payoff;

	// Check all possible next period borrowing choices.
	for (int a_next_index = 0; a_next_index < number_of_a_values; a_next_index++)
	{
		// Get next credit value and the value of consumption implied by the budget constraint.
		double a_next = a_grid[a_next_index];
		double consumption = a_start + e_start - q * a_next;

		// Check budget constraint : if consumption is negative, skip this value.
		if (consumption <= 0.0) { continue; }

		// Calculate the Value Function value.
		current_period_payoff = pow(consumption,(1 - sigma)) / (1 - sigma);

		future_period_payoff = 0;
		for (int e_prob_index = 0; e_prob_index < number_of_e_values; e_prob_index++)
		{
			future_period_payoff += Value_Function[a_next_index][e_prob_index] * e_probs[e_start_index][e_prob_index];
		}
		
		new_v_max = current_period_payoff + beta * future_period_payoff;

		// Check if this capital choice gives the highest Value Function value.
		if (new_v_max > v_max)
		{
			// Update Candidate values.
			v_max = new_v_max;
			a_next_optimal = a_next;
			a_next_optimal_index = a_next_index;
		}
	}

	// Write optimal values.
	V = v_max;
	g = a_next_optimal;
	g_index = a_next_optimal_index;
}

/*
-----------------------------------------------------------------------------------------------------
Function 2 - Solve Value Function and Policy Function.
-----------------------------------------------------------------------------------------------------
Solves for the Value Function and Policy Function using value function iteration.
Applies the Solve Household Problem function to all possible starting states for each iteration.
Search parameters (max iterations, tolerance, etc.) are defined in the function.
Modifies the input initialized arrays since the function cannot return multiple outputs.

Input:
	- Market price for bond
	- Initialized arrays to modify

Output (as changes to input arrays):
	- Value Function
	- Policy Function
	- Policy Function (with indices instead of values)
 */

void Solve_Value_Function(double q, Parameters &params, Array2D<double>& Value_Function,
	Array2D<double>& Policy_Function, Array2D<int>& Policy_Function_Index)
{
	// Unpack relevant parameters.
	vector<double> a_grid = params.a_grid;
	vector<double> e_grid = params.e_grid;

	int number_of_a_values = params.number_of_a_values;
	int number_of_e_values = params.number_of_e_values;

	// Arrays to hold next value function iteration.
	Array2D<double> Value_Function_New(number_of_a_values, vector<double>(number_of_e_values));
	Array2D<double> Policy_Function_New(number_of_a_values, vector<double>(number_of_e_values));
	Array2D<int> Policy_Function_Index_New(number_of_a_values, vector<int>(number_of_e_values));

	// Iteration parameters.
	double dist = DBL_MAX;
	int iteration_count = 0;
	int max_iterations = 5000;
	double tolerance = 1e-6;

	// Optimal values for specific starting state.
	double V, g;
	int g_index;

	// Solve for the Value Function and Policy Function.
	while ((dist > tolerance) && (iteration_count < max_iterations))
	{
		// Loop over all possible starting states.
		for (int a_start_index = 0; a_start_index < number_of_a_values; a_start_index++)
		{
			for (int e_start_index = 0; e_start_index < number_of_e_values; e_start_index++)
			{
				// Solve the Value Function and Policy Function and update values.
				Solve_HH_Problem(q, a_start_index, e_start_index, Value_Function, params, V, g, g_index);

				Value_Function_New[a_start_index][e_start_index] = V;
				Policy_Function_New[a_start_index][e_start_index] = g;
				Policy_Function_Index_New[a_start_index][e_start_index] = g_index;
			}
		}

		// Update search parameters.
		dist = GetDistanceArray2D(Value_Function, Value_Function_New);
		iteration_count++;

		// Update the Value Function and Policy Function.
		DeepCopyArray2D<double>(Value_Function_New, Value_Function);
		DeepCopyArray2D<double>(Policy_Function_New, Policy_Function);
		DeepCopyArray2D<int>(Policy_Function_Index_New, Policy_Function_Index);

		// Print warning if convergence is not achieved.
		if (iteration_count >= max_iterations)
		{
			std::cout << "Warning: value function did not converge after " << iteration_count << " iterations. \n";
		}
	}
}

/*
-----------------------------------------------------------------------------------------------------
Function 3 - Get Population Distribution from Policy Function.
-----------------------------------------------------------------------------------------------------
Solves for the steady-state distribution over credit and endowment states given
a policy function (with index values).
Search parameters (max iterations, tolerance, etc.) are defined in the function.
Modifies an initialized Population Distribution array rather than returning a 
new array to be consistent with Function 1 that returns/modifies multiple arrays.

Input:
	- Policy Function (with index values)
	- Initialized Population Distribution array to modify
Output (as changes to the input array):
	- Steady-state population distribution
 */

void Get_Population_Distribution(Array2D<int>& Policy_Function_Index, Parameters& params,
								 Array2D<double>& Population_Distribution)
{
	// Unpack relevant parameters.
	vector<double> a_grid = params.a_grid;
	vector<double> e_grid = params.e_grid;
	Array2D<double> e_probs = params.e_probs;
	int number_of_a_values = params.number_of_a_values;
	int number_of_e_values = params.number_of_e_values;

	// Array to store next iteration of finding population distribution.
	Array2D<double> New_Distribution(number_of_a_values, vector<double>(number_of_e_values));

	// Initialize with suitable guess distribution.
	double uniform_density = 1.0 / (double(number_of_a_values) * double(number_of_e_values));
	for (int i = 0; i < number_of_a_values; i++)
	{
		for (int j = 0; j < number_of_e_values; j++)
		{
			Population_Distribution[i][j] = uniform_density;
		}
	}

	// Iteration parameters.
	double dist = DBL_MAX;
	int iteration_count = 0;
	int max_iterations = 5000;
	double tolerance = 1e-10;

	// Solve for the steady-state Population Distribution.
	while ((dist > tolerance) && (iteration_count < max_iterations))
	{
		// Get "inflow" to each credit-endowment state in the next period.
		for (int a_index = 0; a_index < number_of_a_values; a_index++)
		{
			for (int e_index = 0; e_index < number_of_e_values; e_index++)
			{
				// Track population inflow into given state.
				double inflow = 0;

				// Loop over all inflow indices.
				for (int a_index_in = 0; a_index_in < number_of_a_values; a_index_in++)
				{
					for (int e_index_in = 0; e_index_in < number_of_e_values; e_index_in++)
					{
						// Check if this state flows into the given state.
						if (Policy_Function_Index[a_index_in][e_index_in] == a_index)
						{
							// Calculate the probability-weighted inflow.
							inflow += Population_Distribution[a_index_in][e_index_in] * e_probs[e_index_in][e_index];
						}
					}
				}
				// Assign total inflow to the Population Distribution.
				New_Distribution[a_index][e_index] = inflow;
			}
		}

		// Update search parameters.
		dist = GetDistanceArray2D(Population_Distribution, New_Distribution);
		iteration_count++;

		// Update the Population Distribution;
		DeepCopyArray2D<double>(New_Distribution, Population_Distribution);

		// Print warning if convergence is not achieved.
		if (iteration_count >= max_iterations)
		{
			std::cout << "Warning: population distribution did not converge after " << iteration_count << " iterations. \n";
		}
	}
}

// Function to calculate the market-clearing condition (sum product of population distribution and policy function arrays).
double Calculate_Market_Clearing_Condition(Array2D<double>& Population_Distribution, Array2D<double>& Policy_Function)
{
	double mcc = 0;
	int n_rows = Population_Distribution.size();
	int n_cols = Population_Distribution[0].size();

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			mcc += Population_Distribution[i][j] * Policy_Function[i][j];
		}
	}
	return mcc;
}

// Function to get the distance between 2 arrays. Defined here as the absolute value of largest element-wise difference.
double GetDistanceArray2D(Array2D<double>& Array1, Array2D<double>& Array2)
{
	double max_diff = -DBL_MAX;
	double element_diff;
	int n_rows = Array1.size();
	int n_cols = Array1[0].size();

	// Loop over all array elements.
	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			// Get the difference between elements.
			element_diff = abs(Array1[i][j] - Array2[i][j]);

			// Check if this is the largest difference.
			if (element_diff > max_diff) { max_diff = element_diff; }
		}
	}
	// Return the largest difference.
	return max_diff;
}

// Function to copy the values in one array to another.
template<class number_type>
void DeepCopyArray2D(Array2D<number_type>& Array_Original, Array2D<number_type>& Array_Copy)
{
	int n_rows = Array_Original.size();
	int n_cols = Array_Original[0].size();

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = 0; j < n_cols; j++)
		{
			Array_Copy[i][j] = Array_Original[i][j];
		}
	}
}

// Function to write a 2-dimensional array to csv file.
template<class number_type>
void WriteArrayToCSV(Array2D<number_type>& Array, const char* file_name)
{
	int n_rows = Array.size();
	int n_cols = Array[0].size();

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

// Instantiate template functions with available types that are not instatiated (called) in this module.
template void WriteArrayToCSV<double>(Array2D<double>& Array, const char* file_name);