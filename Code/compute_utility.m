function discounted_avg_utility = compute_utility(a_val, theta_val, rho_val, num_periods, discount_factor)
     % COMPUTE_UTILITY Computes the discounted average utility by running a Dynare simulation.
    %
    % Parameters:
    %   a_val           - Parameter 'a' value for the simulation.
    %   theta_val       - Parameter 'theta' value for the simulation.
    %   rho_val         - Parameter 'rho' value for the simulation.
    %   num_periods     - Number of periods to consider for utility calculation.
    %   discount_factor - Discount factor for utility.
    %
    % Returns:
    %   discounted_avg_utility - The computed discounted average utility.

    % Declare global variables used by Dynare
    global oo_ M_ options_ estim_params_ bayestopt_

    % Compute discount weights
    weights = discount_factor .^ ((1:num_periods) - 1);

    try
        % Construct Dynare command string with correct formatting.
        % Use 'silent' to reduce output, 'nolog' to avoid logs, 'noclearall' to prevent clearing.
        dynare_command = sprintf('dynare Simulation_loglin noclearall nolog nograph -Da_val=%f -Dtheta_val=%f -Drho_val=%f', a_val, theta_val, rho_val);
        eval(dynare_command);  % Execute the Dynare command

        % Check if Dynare produced the expected IRFs
        if exist('oo_', 'var') && isfield(oo_, 'irfs') && isfield(oo_.irfs, 'U_eps_g')
            % Extract the utility IRF for the specified number of periods
            utility_series = oo_.irfs.U_eps_g(1:num_periods);

            % Compute the discounted average utility
            discounted_avg_utility = sum(weights .* utility_series);
        else
            % Assign a large negative utility if IRFs are not found
            discounted_avg_utility = -1e10;
            warning('Dynare did not produce expected IRFs for a=%.2f, theta=%.2f, rho=%.2f', a_val, theta_val, rho_val);
        end
    catch ME
        % Handle any errors that occur during Dynare execution
        warning('Dynare failed for a=%.2f, theta=%.2f, rho=%.2f: %s', a_val, theta_val, rho_val, ME.message);
        discounted_avg_utility = NaN;
    end
end
