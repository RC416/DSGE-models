/*
Idiosyncratic Endowment model implemented in C++.

Steps:
    1 - Define utility parameters, grids, and parameter struct.
    2 - Solve model and get market bond price using binary search.
        - Guess bond price
        - Solve for the Value Function and Policy Function
        - Get the distribution of credit and productivity levels
        - Check the market clearing condition and update bond price guess
    3 - Plot results.
*/

#include "custom_functions.h"
#include <iostream>
#include <array>

int main()
{
    // -----------------------------------------------------------------------------------------------------
    // 1 - Define utility parameters, grids, and parameter struct.
    // -----------------------------------------------------------------------------------------------------

    // Endowment parameters.
    double e_high = 1.0;                                            // high endowment
    double e_low = 0.1;                                             // low endowment
    int const number_of_e_values = 2;                           
    double e_grid[number_of_e_values] = { e_low, e_high };          // endowment grid
    double e_probs[number_of_e_values * number_of_e_values]
        = { 0.500, 0.500, 0.075, 0.925 };                           // endowment transition probabilities
        
    // Utility parameters.
    double sigma = 1.5;                                             // risk aversion coefficient
    double beta = 0.99322;                                          // discount factor

    // Credit parameters.
    double a_high = 4;                                              // upper credit limit
    double a_low = -2;                                              // lower credit limit
    int const number_of_a_values = 100;                             // credit grid size

    // Create grid of credit values.
    double a_grid[number_of_a_values];
    for (int i = 0; i < number_of_a_values; i++)
    {
        a_grid[i] = a_low + (a_high - a_low) * double(i) / double(number_of_a_values - 1);
    }

    // Store parameters in a struct for passing to a function.
    Parameters params = { sigma, beta, number_of_a_values, number_of_e_values, a_grid, e_grid, e_probs };

    // Value Function, Policy Functions, Population Distribution.
    double Value_Function[number_of_a_values * number_of_e_values] = { 0 };
    double Policy_Function[number_of_a_values * number_of_e_values];
    int Policy_Function_Index[number_of_a_values * number_of_e_values];
    double Population_Distribution[number_of_a_values * number_of_e_values];

    // Placeholders for next Value Function iteration and search for q.
    double Value_Function_New[number_of_a_values * number_of_e_values];
    double Policy_Function_New[number_of_a_values * number_of_e_values];
    int Policy_Function_Index_New[number_of_a_values * number_of_e_values];
    double Population_Distribution_New[number_of_a_values * number_of_e_values];

    // -----------------------------------------------------------------------------------------------------
    // 2 - Solve model and get market bond price using binary search.
    // -----------------------------------------------------------------------------------------------------

    // Range of bond prices values to search.
    double q_min = 0.985;
    double q_max = 1.100;

    // Optional: floor for q_min such that those at credit limit can still afford positive consumption.
    q_min = (a_low + e_low) / a_low;

    // Placeholder for the market clearing condition.
    double mcc = DBL_MAX;

    // Iteration parameters.
    double dist = DBL_MAX;
    int iteration_count = 0;
    int max_iterations = 20;
    double tolerance = 1e-3;

    // Solve for market price q.
    while ((dist > tolerance) && (iteration_count < max_iterations))
    {
        // Get value of q from middle of range.
        double q = (q_min + q_max) / 2;

        // Solve for the Value Function and Policy Function.
        Solve_Value_Function(q, params, Value_Function, Policy_Function, Policy_Function_Index,
            Value_Function_New, Policy_Function_New, Policy_Function_Index_New);

        // Get the Population Distribution.
        Get_Population_Distribution(Policy_Function_Index, params, 
            Population_Distribution, Population_Distribution_New);

        // Check market clearing condition.
        mcc = Calculate_Market_Clearing_Condition(Population_Distribution, Policy_Function,
            number_of_a_values, number_of_e_values);

        // Update search parameters;
        dist = abs(mcc);
        iteration_count++;

        // Update range of q according to the sign of the market clearing condition.
        if (mcc > 0) { q_min = q; };
        if (mcc < 0) { q_max = q; };

        // Print results.
        std::cout << "Iteration " << iteration_count << ": q=" << q << ", mcc=" << mcc << "\n";
        if (iteration_count >= max_iterations)
        {
            std::cout << "Warning: search for q did not converge after " << iteration_count << " iterations \n";
        }

    }

    // -----------------------------------------------------------------------------------------------------
    // 3 - Display results and write to csv files. 
    // -----------------------------------------------------------------------------------------------------

    // Write final Value Function, Policy Function, and Population Distribution to csv files.
    WriteArrayToCSV<double>(Value_Function, "Value_Function.csv", number_of_a_values, number_of_e_values);
    WriteArrayToCSV<double>(Policy_Function, "Policy_Function.csv", number_of_a_values, number_of_e_values);
    WriteArrayToCSV<double>(Population_Distribution, "Population_Distribution.csv", number_of_a_values, number_of_e_values);

    // Display subset of results: Value Function for certain capital and productivity terms.
    std::cout.precision(4);
    std::cout << "\n Value Function \n e = " << e_grid[0] << " \t e = " << e_grid[1] << "\n";
    for (int a_index = 0; a_index < 10; a_index++)
    {
        for (int e_index = 0; e_index < 2; e_index++)
        {
            std::cout << Value_Function[a_index * number_of_e_values + e_index] << "\t \t";
        }
        std::cout << "\n";
    }

	return 0;
}