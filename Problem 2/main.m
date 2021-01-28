%% ASEN 3128 - Lab 1 - Group 4 - Problem 2 - Main
%
% Matlab simulation of the translational dynamics of a ball moving through 
% the air, where the forces on the body are not a function of the body 
% attitude but include aerodynamic drag and gravity.
%
% Authors: R. Charles, C. MacPherson, and Z. Vanlangendonck
% Date: 27th January 2021

%% Housekeeping
clc;
clear;
close all;

%% Constant Declaration

m = 0.03; % mass of golf ball [kg]
d = 0.03; % diameter of golf ball [m]
A = (pi*d^2)/4; % cross-sectional area of golf ball [m^2]
C_d = 0.6; % coefficient of drag
g = 9.81; % acceleration dueto gravity [m/s]
rho = 1.2754; % air density at sea level [Pa]
t_span = [0 5]; % time span for intergration
x_0 = [0 0 0 0 -20 20]; % initial state vector [m/s]
wind = [0 0 0]; % wind vector [m/s]

%% Trajectory Calculations

fprintf("Running trajectory analysis...\n\n");

[t,x] = ode45(@(t,x) objectEOM(t,x,rho,C_d,A,m,g,wind),t_span,x_0);

% Find Landing Location
LL = zeros(1,3);
LLindex = find(x(:,3)<=0,2);
LL(2) = x(LLindex(2),2);
LL(1) = x(LLindex(2),1);

%% Plotting Trajectory

fignum = 1;

fprintf("Plotting trajectory results from ode45 to figure %i...\n\n",fignum);

figure(fignum)
plot3(x(:,1),x(:,2),x(:,3),'linewidth',2);
grid on
title('Golf Ball Trajectory (wind = [0,-20,20][m/s])');
xlabel('x position [m]');
ylabel('y position [m]');
zlabel('z position [m]');
zlim([0 max(x(:,3))*1.1]);
fignum = fignum + 1;

%% Wind Sensitivity Calculations

fprintf("Running wind sensitivity analysis...\n\n");

% initializing vector to store landing locations with changing winds
max_wind = 50;
landing_pos = zeros((max_wind*2)+1,1);

for i = -max_wind:max_wind
    
    % initialize wind vector
    wind = [i 0 0];
    
    [t,x] = ode45(@(t,x) objectEOM(t,x,rho,C_d,A,m,g,wind),t_span,x_0);
    
    index = find(x(:,3)<=0,2);
    landing_pos(i+max_wind+1) = x(index(2),1);
    
end

% finding deflection per m/s of wind
wind_vec = -max_wind:max_wind;
deflect_per_unitWind = zeros((max_wind*2),1);

for i = 1:length(landing_pos)-1
    
    if i <= 50
        deflect_per_unitWind(i) = landing_pos(i)/wind_vec(i);
    else
        deflect_per_unitWind(i) = landing_pos(i+1)/wind_vec(i+1);
    end
    
end

deflect_per_unitWind = mean(deflect_per_unitWind);

fprintf('Found the landing site deflection per m/s of wind --> %0.4f [m/(m/s)] \n\n',deflect_per_unitWind);

%% Plot Wind Sensitivity Findings
 
fprintf("Plotting wind sensittivity results to figure %i...\n\n",fignum);

figure(fignum)
plot(-max_wind:max_wind,landing_pos(:,1),'linewidth',2);
title('Horizontal Wind Sensitivity');
xlabel('Wind Speed [m/s]');
ylabel('Distance from Zero [m]');
fignum = fignum + 1;

%% Kinetic Energy Limitations

fprintf("Running kinetic energy limitation analysis...\n\n");

max_KE = (1/2)*m*(norm(x_0(4:6)))^2; % max kinetic energy [J]
gap = 0.025;
max = 1;
masses = gap:gap:max; % defining varying massed golf ball [kg]

v = sqrt((2*max_KE)./masses)'; % velocities of each mass at given kinetic energy [m/s]

% defining new state vectors for each mass and a vector to record landing
% y-position
x_0_2 = zeros(max/gap,6);
x_0_2(:,5:6) = [-sqrt((v.^2)/2),sqrt((v.^2)/2)];
y_pos = zeros(max/gap,1);

wind = [0 0 0];
t_span = [0 10];

% Find the downrange distances for each mass
for i = 1:(max/gap)
    
    [t,x] = ode45(@(t,x) objectEOM(t,x,rho,C_d,A,masses(i),g,wind),t_span,x_0_2(i,:));
    
    index = find(x(:,3)<=0,2);
    y_pos(i) = abs(x(index(2),2));
    
end

%% Plot Kinetic Energy Limiation Findings

fprintf('Plotting kinetic energy limitation analysis results to figure %i...\n\n',fignum);

figure(fignum)
plot(masses,y_pos,'linewidth',2);
title('Downrange Distance vs Projectile Mass (Constant KE)');
xlabel('Mass of Golf Ball [kg]');
ylabel('Downrange Distance [m]');
fignum = fignum + 1;


