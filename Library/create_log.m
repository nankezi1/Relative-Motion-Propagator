function create_log(log, options, final_distance)
% Description: this function handles the logging of data saving the options
% that were chosen to obtain the related results.

T_DSG = 6.4;                        % days  
tf = options.Nperiods * T_DSG;     % days

fprintf(log, "Settings\n" + ...
             "_________________________________________________________\n\n" + ...
             "Npoints:\t\t%.0f\n" + ...
             "tf:\t\t\t\t%.1f days\n" + ...
             "RelTol:\t\t\t%.13f\n" + ...
             "AbsTol:\t\t\t%.13f\n\n", ...
             options.Npoints, tf, options.RelTolODE, options.AbsTolODE);

fprintf('The Final Distance between Chaser and Target is: %.4f km\n', final_distance);
fprintf(log, '\nThe Final Distance between Chaser and Target is: %.4f km\n', final_distance);


end