%体积
clear

k_vals = 0.05;  % Range of possible wear coefficients
H      = 5e10;  % Hardness [Pa]
Fn     = 686;    % wight [N]
s0     = 0.01;  % Slip distance per footstep [m]
years_used = 500;

W_measured = [0.00344070319919769	0.00319007091714000	0.00248152525640727...
    0.00173072545425254	0.00146184367861409	0.00200781611432773];   % Volume of wear [m^3]

% RESULTS STORAGE
f_steps_array = zeros(size(W_measured));

for i = 1:length(W_measured)
    W_current = W_measured(i);
    % Archardns law => s = V * H / (k * Fn)
    s_total = (W_current * H) / (k_vals * Fn);
    
    % Number of footfalls
    f_steps = s_total / s0;
    f_steps_array(i) = f_steps;
end

% Convert to usage per year (assuming continuous usage over 300 years)
f_steps_per_year = f_steps_array / years_used;

% Display some results
disp('Possible usage rates (footfalls/year) for various wear coefficients:');
disp(table(W_measured', f_steps_per_year', 'VariableNames', ...
           {'W', 'f'}));

%%
%深度
clear

k_vals = 0.05;  % Range of possible wear coefficients
H      = 5e10;  % Hardness [Pa]
Pn     = 30000;    % wight [Pa]
s0     = 0.01;  % Slip distance per footstep [m]
years_used = 500;

W_measured = [2.75494166405037	2.61586063889118	1.75092293274812...
    3.53323352844741	1.25984857145743	2.95858067693070];   % depth of wear [cm]

% RESULTS STORAGE
f_steps_array = zeros(size(W_measured));

for i = 1:length(W_measured)
    W_current = W_measured(i);
    % Archardns law => s = W * H / (k * Pn)
    s_total = (W_current * H) / (k_vals * Pn);
    
    % Number of footfalls
    f_steps = s_total / s0;
    f_steps_array(i) = f_steps;
end

% Convert to usage per year (assuming continuous usage over 300 years)
f_steps_per_year = f_steps_array / years_used;

% Display some results
disp('Possible usage rates (footfalls/year) for various wear coefficients:');
disp(table(W_measured', f_steps_per_year', 'VariableNames', ...
           {'W', 'f'}));
