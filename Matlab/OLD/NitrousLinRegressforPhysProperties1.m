function [ Hv,Cp ] = NitrousLinRegressforPhysProperties1( T )
% NitrousLinRegressProperties
%   The following function uses known values for specific heat capcity with
%   constant pressure, latent heat of vaporization, liquid nitrous density,
%   gaseous nitrous density, and nitrous vapour pressures at a
%   corresponding temperature to find the approximate value at an inputted
%   temperature using linear regression

% Inputs:
% T = Temperature to be evaluated (C)

% Outputs:
% Hv = Latent heat of vaporization at given temp (J/kg)
% Cp = Specific heat capacity w/ const. pressure at given temp (j/kg k)


% corrsponding temperature vector for the given laten heats of vapor,
% specific heat cap., liq. density, vapor density, and vapor pressure. 14
% elements
NitrousTemp = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 36.42]; % Units = C
% Vector for latent heat of vaporization of nitrous oxide. Vector = 14
% elements
HvVal = [285, 276, 266, 255, 244, 232, 219, 204, 188, 169, 147, 117, 64.9, 0]; % Units = Kj/kg
% vector for specific heat capacities, w/ constant pressure process. Last
% two values are made up. At 36.42 (crit Temp) Cp = infinity
CpVal = [1.915,1.957,2.011,2.079,2.166,2.274,2.412,2.592,2.834,3.188,3.781,5.143, 10, 100000]; % Units = Kj/kg K

% Will use linear regression with value of temp
% fall in between two set temperature values, to find the various physical
% properties
% loop control variable to find where the nitrous temperature falls
i = 1;
% ensures that the temperature is within reasonable limits
if T >= 36.42 || T < -25
    disp('Error: Temperature out of bounds');
    display(T);
end
% finds between what given values T falls
while T > NitrousTemp(i)
    i = i+1;
end
% sets Hv equal to HvVal or linear interpolates, does same for CpVal and Cp
if T == NitrousTemp(i)
    Hv = HvVal(i); % units = Kj/kg
    Cp = CpVal(i); % units = Kj/kg K
    
else
    Hv = (((HvVal(i)-HvVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+HvVal(i-1); % units = Kj/kg
    Cp = (((CpVal(i)-CpVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+CpVal(i-1); % units = Kj/kg K
    
end
% Conerts Hv to units of  J/kg
Hv = Hv*1000;
% converts Cp to units of J/kg K
Cp = Cp*1000;

end

