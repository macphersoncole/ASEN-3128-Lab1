%% ASEN 3128 - Lab 1 - objectEOM
% Provides the derivative of the state vector, x, as a function of 
% time, the state, and the rest of the parameters of the problem. 
%
% Authors: R. Charles, C. MacPherson, and Z. Vanlangendonck
% Date: 27th January 2021
function xdot = objectEOM(t,x,rho,C_d,A,m,g,wind)
%
% Inputs:   t    = time [s]
%           x    = state vector 
%                = [x [m], y [m], z [m], x velocity [m/s],
%                   y velocity [m/s], z velocity [m/s]]
%           rho  = density of air [kg/m^3]
%           C_d  = coefficient of drag 
%           A    = cross-sectional area [m^2]
%           m    = mass [kg]
%           g    = acceleration due to gravity [m/s^2]
%           wind = wind velocity vector
%                = [x velocity [m/s], y velocity [m/s], z velocity [m/s]]
%
% Outputs: xdot  = state vector derivatives
%                = [x velocity [m/s], y velocity [m/s], z velocity [m/s], 
%                   x acceleration [m/s], y acceleration [m/s], 
%                   z acceleration [m/s],]
%

%% Initialization

xdot = zeros(6,1);

%% Finding Velocity

% derivative of position is velocity... velocity is given in state equation
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

%% Finding drag

D_x = abs((1/2)*C_d*rho*x(4)*A);
D_y = abs((1/2)*C_d*rho*x(5)*A);
D_z = abs((1/2)*C_d*rho*x(6)*A);

%% Finding Acceleration

xdot(4) = (-D_x/m) + wind(1);
xdot(5) = (-D_y/m) + wind(2);
xdot(6) = (-D_z/m) - g + wind(3);

end