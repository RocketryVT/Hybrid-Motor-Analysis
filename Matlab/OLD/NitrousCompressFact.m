function [ Z ] = NitrousCompressFact( VapPres )
% The following function calculates the compressibility factor of nitrous
% oxide vapor. Assumes the vapor is about saturated, a generally good assumption
% for nitrous emptying out of a tank. Assumes linear line from pressure = 0
% to the critical pressure. Uses linear regression to get value desired

% Inputs
% VapPres = Pressure to find compressibility factor at. % units = bar
% Outputs
% Z = compressiblity factor at given pressure % unitless

% variable for compressibility factor at nitrous oxides critical pressure
Zcrit = 0.28; % unitless
% variable for lower interpolation limit (0) X val;
lowIntpX = 0; % units = bar
% Compressibility factor at pressure = 0 (lowIntpX)
lowIntpY = 1; % unitless
% variable for nitrous oxide critical pressure. Upper Y intp limit
Pcrit = 72.51; % units = bar
% calls interpolation function to calculate Z
Z = LinearInterpolate(VapPres,lowIntpX,lowIntpY,Pcrit,Zcrit); % unitless

end

