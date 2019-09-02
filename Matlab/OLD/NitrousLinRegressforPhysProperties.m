function [ Hv,Cp,LiqD,VapD,VapPres] = NitrousLinRegressforPhysProperties( T )
%NitrousLinRegressProperties Summary of this function goes here
%   The following function uses known values for specific heat capcity with
%   constant pressure, latent heat of vaporization, liquid nitrous density,
%   gaseous nitrous density, and nitrous vapour pressures at a
%   corresponding temperature to find the approximate value at an inputted
%   temperature using linear regression

%Inputs:
    %T = Temperature to be evaluated (C)
    
%Outputs:
    %Hv = Latent heat of vaporization at given temp (J/kg)
    %Cp = Specific heat capacity w/ const. pressure at given temp (j/kg)
    %LiqD = liquid nitrous density at given temp (kg/m^3)
    %VapD = vapor density at given T (kg/m^3)
    %VapPres = Vapor Pressure at given Temp. (bar)
  

%corrsponding temperature vector for the given laten heats of vapor,
%specific heat cap., liq. density, vapor density, and vapor pressure. 14
%elements
NitrousTemp = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 36.42]; %Units = C
%Vector for latent heat of vaporization of nitrous oxide. Vector = 14
%elements
HvVal = [285, 276, 266, 255, 244, 232, 219, 204, 188, 169, 147, 117, 64.9, 0]; % Units = Kj/kg
%vector for specific heat capacities, w/ constant pressure process. Last
%two values are made up. At 36.42 (crit Temp) Cp = infinity
CpVal = [1.915,1.957,2.011,2.079,2.166,2.274,2.412,2.592,2.834,3.188,3.781,5.143, 10, 100000]; %Units = Kj/kg K
%Vector for Nitrous liquid density
LiqDVal = [1014.8,995.4,975.2,953.9,931.4,907.4,881.6,853.5,822.2,786.6,743.9,688,589.4,452]; %units = kg/m^3
%vector for nitrous vapour density
VapDVal = [40.11,46.82,54.47,63.21,73.26,84.86,98.41,114.5,133.9,158.1,190,236.7,330.4,452]; %units = kg/m^3
%vector for vapor pressure (ie tank pressure)
VapPresVal = [1547,1801,2083,2397,2744,3127,3547,4007,4510,5060,5660,6315,7033,7251]; %units = kPa

% Will use linear regression with value of temp
%fall in between two set temperature values, to find the various physical
%properties
%loop control variable to find where the nitrous temperature falls
i = 1;
%ensures that the temperature is within reasonable limits
if T >= 36.42 || T < -25
    disp('Error: Temperature out of bounds');
    display(T);
end
%finds between what given values T falls
while T > NitrousTemp(i)
    i = i+1;
end
%sets Hv equal to HvVal or linear interpolates, does same for CpVal and Cp
if T == NitrousTemp(i)
    Hv = HvVal(i); %units = Kj/kg
    Cp = CpVal(i); %units = Kj/kg K
    LiqD = LiqDVal(i); %units = kg/m^3
    VapD = VapDVal(i); %units = kg/m^3
    VapPres = VapPresVal(i); %units = kPa

else
    Hv = (((HvVal(i)-HvVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+HvVal(i-1); %units = j/kg
    Cp = (((CpVal(i)-CpVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+CpVal(i-1); %units = j/kg K
    LiqD = (((LiqDVal(i)-LiqDVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+LiqDVal(i-1); %units = kg/m^3
    VapD = (((VapDVal(i)-VapDVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+VapDVal(i-1); %units = kg/m^3
    VapPres = (((VapPresVal(i)-VapPresVal(i-1))/(NitrousTemp(i)-NitrousTemp(i-1)))*(T-NitrousTemp(i-1)))+VapPresVal(i-1); %units = kPa
 
end
%Conerts Hv to units of  J/kg
Hv = Hv*1000;
%converts Cp to units of J/kg K
Cp = Cp*1000;
%converts VapPres to units of bar
VapPres= VapPres*0.01; %units = bar

end

