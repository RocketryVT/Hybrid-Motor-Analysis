function [ At,EndD,NozLen,ThroatD ] = NozDimen( mdot,Pc,Mav,Tt,Pa )
% finds the Area of the throat of the conical nozzle to achieve the desired
% combustion chamber pressure based on mass flow rate, average molecular
% exhaust mass, and temperature at throat

% Inputs
% mdot = mass flow rate, % units = kg/sec
% Pc = desired combustion chamber pressure (bar)
% Mav = average Molecular exhaust mass % units = g/mol
% Tt = temperature at nozzle throat % units = K
% Ambient pressure % units = bar
% Outputs
% At = nozzle throat area % units = m^2
% EndD = End Diameter % units = m
% NozLen = length from throat to exit  % units = m
% ThroatD = Throat diameter = % units = m

% value for Runi in J/kmol/k
Runi = 8314.4621; % UNITS = J/kmol/k
% specific heat ratio
k=1.4; % UNITLESS

% finds the component that is multiplied to find At
Atm = sqrt((Runi*Tt)/(Mav*k)); % UNITS = KJ*KELVIN/G
% converts Pc to units of Pa (N/m^2)
Pc = Pc*100000; % units = Pa
% converts Pa to units of Pa
Pa = Pa*100000; % units = Pa
% calculates pressure at nozzle throat
Pt = Pc*((1+((k-1)/2))^(-k/(k-1))); % % UNITS = Pa
% calculates mach number
Nm = sqrt(((2/(k-1))*((((Pc/Pa)^((k-1)/k)))-1))); % units = m/s
% calculates the noxxle throat area need to cause the desired intial
% pressure drop from manifold to combustion chamber
At = (mdot/Pt)*Atm; % units = m^2
% calculates area of exit
Ae = (At/Nm)*(((1+((k-1)/2)*(Nm^2))/((k+1)/2))^((k+1)/(2*(k-1)))); % units = m^2
% calculates EndD
EndD = 2*(sqrt((Ae/pi))); % units = m
% calculates throatD
ThroatD = 2*(sqrt((At/pi))); % units = m
% variable for angle of expansion of nozzle. uses 15 degrees, set conical
% standard, converts to radians
expAngle = (15/180)*pi; % units = rad
% finds the difference in radius of exit and throat
radDif = (EndD/2)-(ThroatD/2); % units = m
% calculates length from throat to exit
NozLen = radDif/(tan(expAngle)); % units = m
end

