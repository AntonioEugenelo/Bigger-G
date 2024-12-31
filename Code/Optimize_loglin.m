close all; clc; clear;

% Set parameters for the search
num_periods = 20;           % Number of periods
discount_factor = 0.997;    % Discount factor

% Parameter grids
theta_values = [6, 8, 9, 10, 11];
rho_values   = [-1000, -1, 0, 1, 1000];

% Store results in a structure
results = struct();

% Suppress optimization display
options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-4);

for t_idx = 1:length(theta_values)
    for r_idx = 1:length(rho_values)
        current_theta = theta_values(t_idx);
        current_rho   = rho_values(r_idx);

        fprintf('----------------------------\n');
        fprintf('Running optimization for theta=%.2f, rho=%.2f\n', current_theta, current_rho);
    
        % Coarse grid over a
        a_values_coarse = linspace(-20, 20, 11);
        Discounted_utilities_coarse = zeros(length(a_values_coarse), 1);

        for i = 1:length(a_values_coarse)
            a_val = a_values_coarse(i);
            Discounted_utilities_coarse(i) = compute_utility(a_val, current_theta, current_rho, num_periods, discount_factor);
        end

        % Find best coarse
        [best_coarse_utility, idx] = max(Discounted_utilities_coarse);
        a_best_coarse = a_values_coarse(idx);

        % Fine grid search around best coarse point
        delta_a = 1.5;  % half-range around best coarse point for fine search
        fine_points = 16; 
        a_values_fine = linspace(a_best_coarse - delta_a, a_best_coarse + delta_a, fine_points);

        Discounted_utilities_fine = zeros(length(a_values_fine), 1);
        for i = 1:length(a_values_fine)
            a_val = a_values_fine(i);
            Discounted_utilities_fine(i) = compute_utility(a_val, current_theta, current_rho, num_periods, discount_factor);
        end

        [best_fine_utility, idx_fine] = max(Discounted_utilities_fine);
        a_best_fine = a_values_fine(idx_fine);

        % Local optimization using fminsearch
        obj_fun = @(a) -compute_utility(a, current_theta, current_rho, num_periods, discount_factor);
        initial_guess = a_best_fine;
        [a_opt, fval_opt] = fminsearch(obj_fun, initial_guess, options);
        final_utility = -fval_opt;

        fprintf('Results for (theta=%.2f, rho=%.2f):\n', current_theta, current_rho);
        fprintf('Coarse Best:  a=%.2f, U=%.5f\n', a_best_coarse, best_coarse_utility);
        fprintf('Fine Best:    a=%.4f, U=%.5f\n', a_best_fine, best_fine_utility);
        fprintf('Local Opt:    a=%.5f, U=%.5f\n', a_opt, final_utility);

        % Store results
        results(t_idx, r_idx).theta = current_theta;
        results(t_idx, r_idx).rho = current_rho;
        results(t_idx, r_idx).a_opt = a_opt;
        results(t_idx, r_idx).utility = final_utility;

        % ---------------- Generate IRF plots for the optimal solution ----------------
        % Re-run the simulation with the optimal a to get IRFs
        global oo_ M_ options_ estim_params_ bayestopt_
        eval(sprintf('dynare Simulation_loglin noclearall nolog -Da_val=%f -Dtheta_val=%f -Drho_val=%f', a_opt, current_theta, current_rho));

        % Find all open figure windows
        figHandles = findall(0, 'Type', 'figure');
        
        % Save each figure with a filename containing the macro variables
        for iF = 1:length(figHandles)
            fH = figHandles(iF);
            figName = sprintf('IRFs_theta_%.2f_rho_%.2f_a_%.5f_fig%d.fig', current_theta, current_rho, a_opt, iF);
            saveas(fH, figName, 'fig');
        end
        
        % Close all figures before next iteration
        close all
    end
end

% After all computations, print a summary of the results
fprintf('\nSummary of Results:\n');
for t_idx = 1:length(theta_values)
    for r_idx = 1:length(rho_values)
        fprintf('theta=%.2f, rho=%.2f: a=%.5f, Utility=%.5f\n', ...
            results(t_idx,r_idx).theta, results(t_idx,r_idx).rho, ...
            results(t_idx,r_idx).a_opt, results(t_idx,r_idx).utility);
    end
end
