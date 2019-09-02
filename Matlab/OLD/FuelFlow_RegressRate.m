function [ mfuelr,r,mdot,portr, BA] = FuelFlow_RegressRate( portr,Cd, Cl, nitrousMdot, a, b,fuelDen,dt)
% The following function calculates fuel mass flow rate and regression

% Inputs
% portr = the current port radius % units = cm
% Cd = combustion chamber diameter % units = cm
% Cl = combustion chamber length % units = cm
% nitrousMdot = nitrous oxide mass flow rate % units = kg/sec
% a = empirical constant for regression rate calculation % unitless
% b^same
% fuelDen = fuel density % units = g/cm^3
% time step % units = sec

% Outputs
% mfuelr = fuel mass flow rat % units =kg/sec
% r = regression rate % units = mm/sec
% mdot = total mass flowrat % units = kg/sec
% portr = port radius. Acts as ref. variable. % units = cm



% converts nitrous mass flow rate to g/s
nitrousMdot = nitrousMdot*1000; % units = g/s

% finds the thickness of the wax layer(cm)
fuelThick = (((Cd./2)-portr)); % UNITS = CM
% check if the chamber is too large
% variable for amount of fuel thickness lost
fLayLost = 0; % units = cm;
% calculates the avg area of a cross section of the port cm^2
portarea = ((((portr+(.5.*fuelThick)).^2).*pi)); % UNITS = CM^2
% calculates the mass flux of the oxidizer (g/sec*cm^2)
G = (nitrousMdot./portarea); % UNITS = G/SEC CM^2
% finds regression rate (mm/sec)
r = (a.*(G.^b)); % UNITS = MM/SEC
% finds how much of the fuel layer is lost
fLayLost = dt*(r/10); % units = cm
% finds the new port radius (increases)
portr = portr+fLayLost; % units = cm
% ensures that the port radius is not bigger than possible
if (2*portr) > Cd
    %  disp('Error: More fuel eroded than possible.');
    % disp(Cd-(portr*2));
end
% burn area cm^2
BA = 2.*pi.*Cl.*portr; % UNITS = CM^2

% fuelmass flow rate g/s found via regression rate
mfuelr = ((r./10).*fuelDen.*BA)/1000; % UNITS = KG/SEC
% calculates mdot
mdot = mfuelr+(nitrousMdot/1000); % units = kg/sec
end

