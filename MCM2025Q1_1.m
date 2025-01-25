% PARAMETERS (Example guesses)
k_vals = 1.0e-3 : 0.5e-3 : 2.0e-3;  % Range of possible wear coefficients
H      = 1.0e9;                    % Hardness [Pa]
weight = 700;                      % Average weight [N] ~ 70 kg * 9.81
foot_area = 0.03;                  % [m^2], approximate foot contact
pn = weight / foot_area;           % Normal pressure [Pa]

% MEASURED OR ESTIMATED
W_measured = 2.5e-4;   % Volume of wear [m^3]
d_step     = 0.001;    % Slip distance per footstep [m]
years_used = 300;

% RESULTS STORAGE
N_steps_array = zeros(size(k_vals));

for i = 1:length(k_vals)
    k_current = k_vals(i);
    % Archard’s law => s = W*H / (k * pn)
    s_total = (W_measured * H) / (k_current * pn);
    
    % Number of footfalls
    N_steps = s_total / d_step;
    N_steps_array(i) = N_steps;
end

% Convert to usage per year (assuming continuous usage over 300 years)
N_steps_per_year = N_steps_array / years_used;

% Display some results
disp('Possible usage rates (footfalls/year) for various wear coefficients:');
disp(table(k_vals', N_steps_per_year', 'VariableNames', ...
           {'WearCoefficient', 'FootfallsPerYear'}));
