% --- Step 1: Define stair geometry ---
width_stair = 1.0;   % [m] total width of the stair
nx = 200;            % Number of discretized points in the x-direction
x_vals = linspace(0, width_stair, nx);

% --- Step 2: Define material & Archard parameters ---
k = 1.0e-3;       % Wear coefficient (dimensionless)
H = 1.0e9;        % Hardness [Pa]

% --- Step 3: Define a Gaussian (or sum of Gaussians) for normal pressure ---
mu1 = 0.5;        % center of Gaussian #1 [m]
sigma1 = 0.05;    % standard deviation [m]
pmax1 = 2.0e5;    % peak pressure [Pa]

% Single Gaussian example:
% p_n = pmax1 * exp( -((x_vals - mu1).^2)/(2*sigma1^2) );

% If two lanes (two Gaussians):
mu2 = 0.2;
sigma2 = 0.05;
pmax2 = 2.0e5;
p_n = pmax1 * exp( -((x_vals - mu1).^2)/(2*sigma1^2) ) ...
    + pmax2 * exp( -((x_vals - mu2).^2)/(2*sigma2^2) );

% --- Step 4: Relate sliding distance s(x) to traffic ---
% For simplicity, assume s(x) is proportional to p_n(x)
% or define s(x) separately if you have a more detailed model.
C_slip = 1.0;  % some proportionality constant
s_x = C_slip * p_n;  % Example: local sliding distance scales with local traffic

% --- Step 5: Compute local wear w(x) and integrate for total worn volume ---
w_x = k * (p_n .* s_x) / H;  % Archard’s Law in discrete form
dx = x_vals(2) - x_vals(1);
W_total = sum(w_x) * dx;     % approximate integral over x

% --- Step 6: Visualize or analyze the lateral wear profile ---
figure;
subplot(2,1,1);
plot(x_vals, p_n, 'r', 'LineWidth', 2);
xlabel('Position x across stair [m]');
ylabel('Normal Pressure p_n [Pa]');
title('Gaussian Distribution of Pressure');

subplot(2,1,2);
plot(x_vals, w_x, 'b', 'LineWidth', 2);
xlabel('Position x across stair [m]');
ylabel('Local Wear Rate w(x)');
title('Wear Profile from Archard + Gaussian Model');

% Display total worn volume
fprintf('Total worn volume (per unit step depth) = %.3e m^3\n', W_total);
