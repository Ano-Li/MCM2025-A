%体积
clear

k_vals = 0.05;  % Range of possible wear coefficients
H      = 5e10;  % Hardness [Pa]
Fn     = 686;    % wight [N]
s0     = 0.01;  % Slip distance per footstep [m]
years_used = 500;

V_measured = [0.00344070319919769	0.00319007091714000	0.00248152525640727...
    0.00173072545425254	0.00146184367861409	0.00200781611432773...
    0.00102799793499907	0.000789410771073454	0.000787570697464598	0.00142666029777259	0.000914022931391297];   % Volume of wear [m^3]

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
           {'V', 'f'}));


%%
%画图
plotmatrix = [V_measured;f_steps_per_year];
pm_sorted = sort(V_measplotmatrixured); % 从小到大排序
