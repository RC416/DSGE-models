module custom_functions
    implicit none

    ! Data structure to store parameters.
    type Parameters
        real(8) :: alpha, beta, delta
        integer :: number_of_k_values, number_of_z_values
        real(8), allocatable, dimension(:)      :: k_values
        real(8), allocatable, dimension(:)      :: z_values
        real(8), allocatable, dimension(:,:)    :: z_probs
    end type Parameters

    ! Function declarations.
    contains

    ! Function to solve the household's problem for a given starting state.
    subroutine Solve_HH_Problem(Value_Function, Policy_Function, iteration, kt0_index, zt_index, params)

        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.

        ! Input variables.
        real(8), dimension(:,:,:), intent(inout)  :: Value_Function, Policy_Function
        integer, intent(in) ::  iteration, kt0_index, zt_index
        type(Parameters), intent(in) :: params

        ! Variables for solving the household's problem.
        real(8) :: v_max, kt1_optimal, new_value_function_value, kt0, zt
        real(8) :: new_v_max, kt1, current_period_payoff, future_period_payoff
        integer :: kt1_index, z_prob_index

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Get starting capital and productivity values from index.
        kt0 = params%k_values(kt0_index)
        zt = params%z_values(zt_index)

        ! Variables to store candidate optimal values.
        v_max = -huge(kt0)
        kt1_optimal = 0.0
        new_value_function_value = 0.0

        do kt1_index = 1, params%number_of_k_values
            
            ! Get capital value from index.
            kt1 = params%k_values(kt1_index)

            ! Calculate value function for a given next period capital choice.
            current_period_payoff = log(zt * (kt0 ** params%alpha) + (1 - params%delta) * kt0 - kt1)

            future_period_payoff = 0.0

            do z_prob_index = 1, params%number_of_z_values
                future_period_payoff = future_period_payoff + &
                    Value_Function(z_prob_index, kt1_index, iteration - 1) * params%z_probs(zt_index, z_prob_index)
            end do

            new_v_max = current_period_payoff + params%beta * future_period_payoff

            ! Check if this capital choice gives the highest Value Function value.
            if (new_v_max > v_max) then
                
                ! Update candidate values.
                v_max = new_v_max
                kt1_optimal = kt1_index
            end if
        end do

        ! Update the Value Function and Policy Funtion with optimal values.
        Value_Function(zt_index, kt0_index, iteration) = v_max
        Policy_Function(zt_index, kt0_index, iteration) = kt1_optimal
    end subroutine


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


    ! Function to read a vector from csv file (stored column-wise).
    subroutine ReadVectorFromCSV(file_path, n_rows, Output_Vector)
        
        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.
        
        ! Input and output variables.
        character(len=*), intent(in)            ::  file_path
        integer, intent(in)                     ::  n_rows
        real(8), dimension(n_rows), intent(out) ::  Output_Vector
        
        ! Local variables.
        integer ::  i

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Open file.
        open(unit=10, file=file_path, access="sequential")

        ! Read file row-by-row.
        do i = 1, n_rows
            read(10,*) Output_Vector(i)
        end do
        
        ! Close file.
        close(unit=10)

    end subroutine


    ! Function to read a 2-dimensional array from csv file.
    subroutine ReadArrayFromCSV(file_path, n_rows, n_cols, Output_Array)
        
        ! -----------------------------------------------------------------------------------------------
        ! Declaration section.
        
        ! Input and output variables.
        character(len=*), intent(in)                    ::  file_path
        integer, intent(in)                             ::  n_rows, n_cols
        real(8), dimension(n_rows, n_cols), intent(out) ::  Output_Array
        
        ! Local variables.
        integer ::  i

        ! -----------------------------------------------------------------------------------------------
        ! Execution section.

        ! Open file.
        open(unit=10, file=file_path, access="sequential")

        ! Read file row-by-row.
        do i = 1, n_rows
            read(10,*) Output_Array(i,:)
        end do
        
        ! Close file.
        close(unit=10)

    end subroutine

end module custom_functions