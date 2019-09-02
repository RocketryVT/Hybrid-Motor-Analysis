function [ Pc ] = ChamberPressure( mdot,At,Tt,Mav, k )
% Following Function calculates combustion chamber pressure

% Inputs
% nitrousMdot = nitrous mass flow rate % units = kg/sec
% mfuelr = mass flow rate of fuel % units = kg/sec
% Tt = Temperature at nozzle throat % units = K
% Mav = Average Molecular mass of exhaust gases % units = g/mol

% Outputs
% Pc = combustion chamber pressure % units = bar

% value for Runi in J/kmol/k
Runi = 8314.4621; % UNITS = J/kmol/k
% specific heat ratio
% UNITLESS

% finds the component that is multiplied to find At
Atm = sqrt((Runi*Tt)/(Mav*k)); % UNITS = KJ*KELVIN/G
% calculates the combustion chamber pressure based on the nozzle throat area
Pc = ((((mdot)*Atm)/At)/((1+((k-1)/2))^(-k/(k-1)))); % units = Pa
Cstar = 1700; % units = m/s
% Pc = (Cstar*mdot/At);
% converts Pc to kPa
Pc = Pc/1000; % units = kPa

% converts kPa to bar
Pc=Pc*0.01; % units = bar;


end

