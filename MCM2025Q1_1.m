clear

k_vals = 0.03;  % Range of possible wear coefficients
H      = 5e10;  % Hardness [Pa]
Fn     = 686;    % wight [N]
s0     = 0.01;  % Slip distance per footstep [m]
years_used = 500;

V_measured = [0.220205004748652	0.204164538696960	0.158817616410066...
    0.110766429072162	0.0935579954313017	0.128500231316975];   % Volume of wear [m^3]

% RESULTS STORAGE
f_steps_array = zeros(size(V_measured));

for i = 1:length(V_measured)
    V_current = V_measured(i);
    % Archardns law => s = V * H / (k * Fn)
    s_total = (V_current * H) / (k_vals * Fn);
    
    % Number of footfalls
    f_steps = s_total / s0;
    f_steps_array(i) = f_steps;
end

% Convert to usage per year (assuming continuous usage over 300 years)
f_steps_per_year = f_steps_array / years_used;

% Display some results
disp('Possible usage rates (footfalls/year) for various wear coefficients:');
disp(table(V_measured', f_steps_per_year', 'VariableNames', ...
           {'W', 'f'}));
