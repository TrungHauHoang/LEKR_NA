clear; clc; close all;

%% Parameters
L = 1;                % Domain length
Nx = 100;             % Number of grid points
h = L / (Nx);     % Grid spacing
x = linspace(0, L, Nx);
T = 0.2;              % Final time
dt = 0.0001;           % Time step
Nt = round(T/dt);     % Number of time steps

v = 1;                % Advection speed
kappa_physical = 1/(1000);   % Physical diffusion coefficient
kappa_fixed = 0.01;
C = 0.5;              % Artificial diffusion coefficient

%% Initial Condition: Sharp Pulse
U = exp(-200 * (x - 0.5).^2);  
% U = 4*x.*(1-x);  
peclet_ori = (v*h)/(kappa_physical)
peclet_mod = (v*h)/(kappa_physical + kappa_fixed)
CFL_ori = (kappa_physical*dt)/(h^2)
CFL_mod = ((kappa_physical + kappa_fixed)*dt)/(h^2)
%% Apply homogeneous Dirichlet BC: U(0) = U(L) = 0
U(1) = 0;
U(end) = 0;

%% Solve with different diffusion approaches
U_no_diff = U;        % No artificial diffusion
U_fixed_kappa = U;    % Increased kappa from the start
U_adaptive = U;       % Adaptive artificial diffusion

%% Time-stepping loop
for n = 1:Nt
    % Compute numerical derivatives (excluding boundaries)
    
    Ux = (U_no_diff(3:end) - U_no_diff(1:end-2)) / (2*h);  % Centered difference for Ux
    Uxx = (U_no_diff(3:end) - 2 * U_no_diff(2:end-1) + U_no_diff(1:end-2)) / (h^2);  % Centered difference for Uxx
    
    % No artificial diffusion (pure advection)
    U_no_diff(2:end-1) = U_no_diff(2:end-1) - dt * v * Ux + dt * kappa_physical * Uxx;
    
    % Fixed increased kappa
      % Fixed diffusion coefficient
    U_fixed_kappa(2:end-1) = U_fixed_kappa(2:end-1) - dt * v * Ux + dt *(kappa_physical + kappa_fixed) * Uxx;
    
    % Adaptive artificial diffusion
    U_max = max(abs(U_adaptive));  % Approximate max U for artificial diffusion scaling
    kappa_artificial = C * U_max * h / 2;  % Artificial diffusion
    U_adaptive(2:end-1) = U_adaptive(2:end-1) - dt * v * Ux + dt * (kappa_physical + kappa_artificial) * Uxx;
    
    % Enforce homogeneous Dirichlet BCs at each step
    U_no_diff(1) = 0; U_no_diff(end) = 0;
    U_fixed_kappa(1) = 0; U_fixed_kappa(end) = 0;
    U_adaptive(1) = 0; U_adaptive(end) = 0;
end

%% Plot results
figure;
plot(x, U, 'k--', 'LineWidth', 1.5); hold on;
plot(x, U_no_diff, 'r', 'LineWidth', 1.5);
% plot(x, U_fixed_kappa, 'b', 'LineWidth', 1.5);
% plot(x, U_adaptive, 'g', 'LineWidth', 1.5);
legend('Initial', 'No Artificial Diffusion', 'Fixed Increased \kappa', 'Adaptive Artificial Diffusion');
xlabel('x'); ylabel('U');
title('Advection-Diffusion with Artificial Diffusion & Homogeneous Dirichlet BC');
grid on;
