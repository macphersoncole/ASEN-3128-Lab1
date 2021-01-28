% Function for ode45 to utilize while computing
%
% Inputs: t = time
%         pos = position vector
%
% Outputs dpos = end position vector for one iteration of ode45
%

function dydt = odefun(t,y)
    
    dydt = zeros(3,1);
    dydt(1) = y(1) + 2*y(2) + y(3);
    dydt(2) = y(1) - 5*y(3);
    dydt(3) = y(1)*y(2) - y(2)^2 + 3*(y(3)^3);
    
end