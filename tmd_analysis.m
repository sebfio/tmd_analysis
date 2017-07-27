% Taipei 101 TMD Spring Damper Analysis
%
% Sebastian Fiorini
% ID: 20558540
% Date: July 22nd 2017
clc;
clear all;
close all;

L = 44.56;
m = 660E3;
g = 9.81;

k = 111.88E6;
C = 2.0E6;

starting_offset = 1.50 / 2;
starting_ang_velocity = 0;

x1 = 1;
x2 = C / m;
x3 = k / m;

% Homogenous Solution

% Change of variables
% use y1 and y2
syms y(t)

[V] = odeToVectorField(diff(y, 2) + x2 * diff(y) + x3 * y== 0);
M = matlabFunction(V,'vars', {'t','Y'});
sol = ode45(M,[0 5],[starting_offset starting_ang_velocity]);

disp_sol = @(x)deval(sol,x,1);

fplot(disp_sol, [0, 5]);

hold on;
title('Displacement plot of 75cm Offset Starting Condition');
xlabel('Time (s)');
ylabel('Displacment (m)');
grid on;
grid minor;
hold off;

figure;
vel_sol = @(x)deval(sol,x,2);
fplot(vel_sol, [0, 5]);

hold on;
title('Velocity plot of 75cm Offset Starting Condition');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;
grid minor;
hold off;

% Obtain the force applied at the spring damper by multiplying the numerical
% solution for displacement at L\theta by k and L\dot{theta} by C
force_applied = @(x)(C * L * vel_sol(x) + k * L * disp_sol(x));

figure;
fplot(force_applied, [0, 5]);

hold on;
title('Force Response to 150cm Offset Starting Condition');
xlabel('Time (s)');
ylabel('Force (N)');
grid on;
grid minor;
hold off;

% Total Solution - including steady state 
z_0 = -0.75;
Y   = 0.75;
v_0 = 0;

w_n = sqrt(k / m);
zeta = C / (2 * w_n * m);

ohm = w_n;

r = ohm / w_n;
psi = pi / 2;
Z = starting_offset / 2 / zeta / r;
w_d = w_n * sqrt(1 - zeta^2);

phi_d = atan((v_0 + (z_0 - Z * cos(psi)) * zeta * w_n - Z * w_n) / (w_d * (z_0 - Z * cos(psi))));

Z_0 = (z_0 - Z * cos(psi)) / cos(phi_d);

% Deriving components of TMD response
z_steady = @(t) Z * cos(w_n * t - psi);
z_homogenous = @(t) Z_0 * exp(-zeta * w_n * t) * cos(w_d * t - phi_d);

z_full = @(t)(z_homogenous(t) + z_steady(t));
y_full = @(t) Y * cos(13.02 * t);

x_steady    = @(t) y_full(t) + z_steady(t); 
x_full      = @(t) z_full(t) + y_full(t);

% Plot 2 periods of Transient response
plot_times.start    = 0;
plot_times.end      = 2 * (2 * pi) / w_n;

figure;
fplot(x_full, [plot_times.start, plot_times.end]);

hold on;
title('2 Periods of Displacement Response to Base Excitation');
xlabel('Time (s)');
ylabel('Displacement (m)');
ylim([-3 3]);
grid on;
grid minor;
saveas(gcf, 'generated_images/2_periods_base_excitation.jpg');
hold off;

% Plot 1 period of Steady State response
plot_times.end = 1 * (2 * pi) / w_n;

figure;
fplot(x_steady, [plot_times.start, plot_times.end]);

hold on;
title('1 Period Steady State Response for TMD');
xlabel('Time (s)');
ylabel('Displacement (m)');
ylim([-4 4]);
grid on;
grid minor;
saveas(gcf, 'generated_images/1_period_steady_state.jpg');
hold off;


% Find force from DAMPERS during steady state response.
plot_times.end = 1 * (2 * pi) / w_n;

figure;

% NOTE: Function handles can't be easily derivated in matlab, only symbolics
force_damper = @(t) -C * Z * w_n * sin(w_n * t - psi);
force_spring = @(t) k * z_steady(t);

fplot(force_damper, [plot_times.start, plot_times.end]);
hold on;
title('Force Applied Through Damper of TMD');
xlabel('Time (s)');
ylabel('Force (N)');
grid on;
grid minor;
saveas(gcf, 'generated_images/force_from_damper.jpg');
hold off;


figure;
fplot(force_spring, [plot_times.start, plot_times.end]);

hold on;
title('Force Applied Through Spring of TMD');
xlabel('Time (s)');
ylabel('Force (N)');
grid on;
grid minor;
saveas(gcf, 'generated_images/force_from_spring.jpg');
hold off;

figure;
fplot(@(t) force_spring(t) + force_damper(t), [plot_times.start, plot_times.end]);

hold on;
title('Total Force Applied on TMD');
xlabel('Time (s)');
ylabel('Force (N)');
grid on;
grid minor;
saveas(gcf, 'generated_images/total_force_TMD.jpg');
hold off;
