module custom_functions
    implicit none

    ! Data structure to store parameters.
    type Parameters
        real(8) :: sigma, beta
        integer :: number_of_a_values, number_of_e_values
        real(8), allocatable, dimension(:)      :: a_grid
        real(8), allocatable, dimension(:)      :: e_grid
        real(8), allocatable, dimension(:,:)    :: e_probs
    end type Parameters

    ! Function declarations.
    contains

    ! -----------------------------------------------------------------------------------------------------
    ! Function 1 - Solve Household Problem.
    ! -----------------------------------------------------------------------------------------------------
    ! Function to find the next iteration of the Value Function and Policy Function
    ! by solving the household's problem for a given starting state.
    ! Outputting multiple values is handled by passing references to output variables.
    !
    ! Input:
    !   - Current Value Function
    !   - Market price for bond
    !   - Starting wealth and endowment level
    !   - References to output variables (V, g, g_index)
    !
    ! Output:
    !   - Next value function value
    !   - Optimal choice of savings/borrowing for next period

    subroutine Solve_HH_Problem(q, a_start_index, e_start_index, params, Value_Function, V, g, g_index)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.

        ! Input and output variables.
        real(8), intent(in) :: q
        integer, intent(in) :: a_start_index, e_start_index
        type(Parameters), intent(in) :: params  
        real(8), dimension(:,:), intent(in) :: Value_Function
        real(8), intent(inout) :: V, g
        integer, intent(inout) :: g_index

        ! Variables for solving the household's problem.
        real(8) :: v_max, new_v_max, a_next_optimal, a_start, e_start, a_next 
        real(8) :: consumption, current_period_payoff, future_period_payoff
        integer :: a_next_index, a_next_optimal_index, e_prob_index

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Get starting capital and productivity values from index.
        a_start = params%a_grid(a_start_index)
        e_start = params%e_grid(e_start_index)

        ! Variables to store candidate optimal values.
        v_max = -huge(a_start)
        a_next_optimal = -2.0
        a_next_optimal_index = 1

        do a_next_index = 1, params%number_of_a_values
            
            ! Get capital value from index.
            a_next = params%a_grid(a_next_index)

            ! Get the value of consumption implied by the budget constraint.
            consumption = a_start + e_start - q * a_next

            ! Check budget constraint : if consumption is negative, skip this value.
            if (consumption <= 0.0) cycle

            ! Calculate the Value Function value.
            current_period_payoff = (consumption ** (1 - params%sigma)) / (1 - params%sigma)

            future_period_payoff = 0.0
            do e_prob_index = 1, params%number_of_e_values
                future_period_payoff = future_period_payoff + &
                Value_Function(e_prob_index, a_next_index) * params%e_probs(e_start_index, e_prob_index)
            end do

            new_v_max = current_period_payoff + params%beta * future_period_payoff

            ! Check if this capital choice gives the highest Value Function value.
            if (new_v_max > v_max) then
                ! Update candidate values.
                v_max = new_v_max
                a_next_optimal = a_next
                a_next_optimal_index = a_next_index
            end if

        end do
        
        ! Write optimal values
        V = v_max
        g = a_next_optimal
        g_index = a_next_optimal_index

    end subroutine

    ! -----------------------------------------------------------------------------------------------------
    ! Function 2 - Solve for the Value Function and Policy Function.
    ! -----------------------------------------------------------------------------------------------------
    ! Solves for the Value Function and Policy Function using value function iteration.
    ! Applies the Solve Household Problem function to all possible starting states for each iteration.
    ! Search parameters (max iterations, tolerance, etc.) are defined in the function.
    !
    ! Input:
    !    - Market price for bond
    !    - Initialized arrays to modify
    !
    ! Output (as changes to input arrays):
    !    - Value Function
    !    - Policy Function
    !    - Policy Function (with indices instead of values)

    subroutine Solve_Value_Function(q, params, Value_Function, Policy_Function, Policy_Function_Index, &
    Value_Function_New, Policy_Function_New, Policy_Function_Index_New)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.
        
        ! Input and output variables.
        real(8), intent(in) :: q
        type(Parameters), intent(in) :: params
        real(8), dimension(:,:), intent(inout) :: Value_Function, Policy_Function, Value_Function_New, Policy_Function_New
        integer, dimension(:,:), intent(inout) :: Policy_Function_Index, Policy_Function_Index_New

        ! Iteration parameters.
        real(8) :: dist, tolerance
        integer :: iteration_count, max_iterations, a_start_index, e_start_index, i, j

        ! Optimal values for specific starting state.
        real(8) :: V, g
        integer :: g_index

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Set first iteration of Value Function to 0.
        do i = 1, params%number_of_e_values
            do j = 1, params%number_of_a_values
                Value_Function(i,j) = 0
            end do
        end do

        ! Iteration parameters.
        dist = huge(dist);
        iteration_count = 0;
        max_iterations = 5000;
        tolerance = 1e-6;

        ! Solve for the Value Function and Policy Function.
        do while ((dist > tolerance) .and. (iteration_count < max_iterations))

            ! Loop over all possible starting states.
            do a_start_index = 1, params%number_of_a_values
                do e_start_index = 1, params%number_of_e_values

                    ! Solve Value Function and Policy Function and update values.
                    call Solve_HH_Problem(q, a_start_index, e_start_index, params, Value_Function, V, g, g_index)

                    Value_Function_New(e_start_index, a_start_index) = V
                    Policy_Function_New(e_start_index, a_start_index) = g
                    Policy_Function_Index_New(e_start_index, a_start_index) = g_index
                end do
            end do

            ! Update search parameters.
            dist = GetDistanceArray2D(Value_Function, Value_Function_New)
            iteration_count = iteration_count + 1

            ! Update the Value Function and Policy Function.
            Value_Function = Value_Function_New
            Policy_Function = Policy_Function_New
            Policy_Function_Index = Policy_Function_Index_New

            ! Print warning if convergence is not achieved.
            if (iteration_count >= max_iterations) then
                write(*,*) "Warning: value function did not converge after", iteration_count, &
                 "iterations", NEW_LINE('A')
            end if
        end do


    end subroutine

    ! -----------------------------------------------------------------------------------------------------
    ! Function 3 - Get Population Distribution from Policy Function.
    ! -----------------------------------------------------------------------------------------------------
    ! Solves for the steady state distribution over credit and endowment states given
    ! a policy function (with index values).
    ! Search parameters (max iterations, tolerance, etc.) are defined in the function.
    ! Modifies an initialized Population Distribution array rather than returning a 
    ! new array to be consistent with Function 1 that returns/modifies multiple arrays.
    !
    ! Input:
    !    - Policy Function (with index values)
    !    - Initialized Population Distribution array to modify
    ! Output (as changes to the input array):
    !    - Steady-state population distribution

    subroutine Get_Population_Distribution(Policy_Function_Index, params, & 
    Population_Distribution, Population_Distribution_New)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.
        
        ! Input and output variables.
        integer, dimension(:,:), intent(in) :: Policy_Function_Index
        type(Parameters), intent(in) :: params
        real(8), dimension(:,:), intent(inout) :: Population_Distribution, Population_Distribution_New

        ! Iteration variables.
        real(8) :: uniform_density, dist, tolerance, inflow
        integer :: iteration_count, max_iterations, a_index, a_index_in, e_index, e_index_in, i, j

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Initialize with suitable guess distribution.
        uniform_density = 1.0 / (params%number_of_a_values * params%number_of_e_values)
        do i = 1, params%number_of_a_values
            do j = 1, params%number_of_e_values
                Population_Distribution(j, i) = uniform_density
            end do
        end do

        ! Iteration parameters.
        dist = huge(dist);
        iteration_count = 0;
        max_iterations = 5000;
        tolerance = 1e-8;

        ! Solve for the steady-state Population Distribution.
        do while ((dist > tolerance) .and. (iteration_count < max_iterations))

            ! Get "inflow" to each credit-endowment state in the next period.
            do a_index = 1, params%number_of_a_values
                do e_index = 1, params%number_of_e_values
                    
                    ! Track population inflow into given state.
                    inflow = 0

                    ! Loop over all inflow indices.
                    do a_index_in = 1, params%number_of_a_values
                        do e_index_in = 1, params%number_of_e_values

                            ! Check if this state flows into the given state.
                            if (Policy_Function_Index(e_index_in, a_index_in) == a_index) then

                                ! Calculate the probability-weighted inflow.
                                inflow = inflow + Population_Distribution(e_index_in, a_index_in) * &
                                params%e_probs(e_index_in, e_index)
                            end if
                        end do
                    end do

                    ! Assign total inflow to the Population Distribution.
                    Population_Distribution_New(e_index, a_index) = inflow
                end do
            end do

            ! Update search parameters.
            dist = GetDistanceArray2D(Population_Distribution, Population_Distribution_New)
            iteration_count = iteration_count + 1

            ! Update the Population Distribution.
            Population_Distribution = Population_Distribution_New
            
            ! Print warning if convergence is not achieved.
            if (iteration_count >= max_iterations) then
                write(*,*) "Warning: population distribution did not converge after", iteration_count, & 
                 "iterations", NEW_LINE('A')
            end if
        end do

    end subroutine

    !  Function to calculate the market-clearing condition (sum product of population distribution and policy function arrays).
    subroutine Calculate_Market_Clearing_Condition(mcc, Population_Distribution, Policy_Function)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.
        
        real(8), dimension(:,:), intent(in) :: Population_Distribution, Policy_Function
        real(8), intent(inout) :: mcc
        integer :: i, j

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        mcc = 0

        do i = 1, size(Population_Distribution, dim=1)
            do j = 1, size(Population_Distribution, dim=2)
                mcc = mcc + Population_Distribution(i, j) * Policy_Function(i, j)
            end do
        end do
    end subroutine

    ! Function to get the distance between 2 arrays. Defined here as the absolute value of largest element-wise difference.
    function GetDistanceArray2D(Array1, Array2) result(max_diff)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.

        real(8), dimension(:,:), intent(in) :: Array1, Array2
        real(8) :: max_diff, element_diff
        integer :: i, j

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Initialize maximum difference to low value.
        max_diff = -huge(max_diff)

        do i = 1, size(Array1, dim=1)
            do j = 1, size(Array1, dim=2)
                
                ! Get the difference between elements.
                element_diff = abs(Array1(i,j) - Array2(i,j))

                ! Check if this is the largest difference.
                if (element_diff > max_diff) max_diff = element_diff
            end do
        end do
    end function

    ! Function to write a 2-dimensional array to csv file.
    subroutine WriteArrayToCSV(Array2D, n_rows, n_cols, file_name)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.

        ! Input variables.
        integer, intent(in)                             :: n_rows, n_cols
        character(len=*), intent(in)                    :: file_name
        real(8), dimension(n_rows, n_cols), intent(in)  :: Array2D

        ! Local variables. 
        integer :: i

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.
        
        ! Formatting for CSV (label = 101)
        101 format(1x, *(g0, ", "))

        ! Create and open file.
        open(unit = 10, access = "sequential", action = "write", &
            status = "replace", file = file_name, form = "formatted")

        ! Write each row to file.
        do i = 1, n_rows
            write(10, 101) Array2D(i,:)
        end do

        ! Close file.
        close(10)
    
    end subroutine

end module custom_functions