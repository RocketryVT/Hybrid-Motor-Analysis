function [ T,Mv,Mold,Mtot,Mnew,OMv,MvLag,nitrousMdot,MassVap,VapPres,LiqD,VapD] = LiqNitrousEmpty( T,Mv,Mold,TankV,Pc,N,injD,dt,Mtot,Mnew,OMv,MvLag,nitrousMdot,MassVap )
% The following Function is used to simulate liquid nitrous oxide flowing
% out of a tank. It returns so many variable because that way it acts as
% pass by references, thus allowing this function to be used in a loop and
% pick up where it left off
% Developed using Aspire Space Technical paper detailing liquid nitrous
% emptying numerical models
% Inputs
% T = Temp of nitrous (C)
% Mv = Mass of nitrous oxide vaporised in given time step % units = kg
% Mold = Previous mass of liq. nitrous not accounting for vaporization
% % units =kg
% TankV = Volume of nitrous tank % units = m^3
% Pc = combustion chamber pressure % units = bar
% N = number of injector orifices % unitless
% injD = injector orifice diameter % units = mm
% dt = time step % units = sec
% Mtot = total nitrous mass in tank % units = kg
% Mnew = liq nitrous mass after accounting for vaporization % units = kg
% OMv = old mass of nitrous vaporized % units = kg
% MvLag = lagged mass of vaporized nitrous % units =kg
% MassVap = mass of nitrous vapor 5units = kg

% Outputs
% same as input, but changed value, like a pass by reference
% VapPres = Vapor Pressure % units = bar

% variable for loss coeffcient. Empirical. Default =2
K = 2; % unitless
% Pc = 1;
% O means old
% variable for diameter of injector
% injD = 10; % units = mm
% calculates the area of each injector
injA = (((injD/1000)*0.5)^2)*pi; % units = m^2
% expression for the denominator on the mass flow rate eq from aspire space
Dloss = (K/((N*injA)^2))*(1/100000); % units = 1/m^4
% sets the a new old nitrous mass flow rate based on the current nitrousmdot
OnitrousMdot = nitrousMdot; % units = kg/s
% updates nitrous physical properties based on temp
% for function called below
% Outputs:
% Hv = Latent heat of vaporization at given temp (J/kg)
% Cp = Specific heat capacity w/ const. pressure at given temp (J/kg k)

[ Hv,Cp] = NitrousLinRegressforPhysProperties1( T ); % units above
% calculates the heat removed from the liquid nitrous
% OMv =  Old Mass of nitrous vaporized last. NOT Mass of vapor
dQ = OMv*Hv; % Units = J
% finds the change in temperature Mnew is mass of liq. nitrous after
% accounting for vaporization

dT = -(dQ/(Mnew*Cp)); % units = K or C (increment is the same)
% calculates the new nitrous temp
T = T+dT; % units = C

% updates nitrous physcal properties, funt outputs below
% Outputs:

% LiqD = liquid nitrous density at given temp (kg/m^3)
% VapD = vapor density at given T (kg/m^3)
% VapPres = Vapor Pressure at given Temp. (bar)
[LiqD,VapD,VapPres ] = NitrousLinRegressforPhysProperties2( T ); % units above

% calculates pressure difference between chamber and Nitrous
dPres = VapPres-Pc; % units = bar
% dPres=10;
% calculates nitrous mdot
nitrousMdot = sqrt(((2*LiqD*dPres)/Dloss)); % units = kg/sec
% calculates the mass lost via Addams Second Order Integration forumla:
% Xn = X(n-1)+dt/2 * ((3*Xdot(n-1)-Xdot(n-2))
MassLost = 0.5*dt*(3*nitrousMdot-OnitrousMdot); % units = kg
% calculate the new total mass in the tank
Mtot = Mtot-MassLost; % units = kg
% the following variable is the mass of liquid nitrous prior to vaporization
% effects
Mold = Mold-MassLost; % units = kg
% Now the Effects of vaporization are accounted for-------------
% calculates the mass of liquid nitrous accounting for vaporization
Mnew = (TankV-((Mtot/VapD)))/((1/LiqD)-(1/VapD)); % units = kg
% calculates the mass of nitrous vapor
MassVap = Mtot-Mnew; % units = kg
% Now calculates the amount of nitrous vaporized in the last itme interval
% Mv = mass vaporized (v = vaporized)
Mv = Mold - Mnew; % units = kg
% the following adds a first order time lag (0.15 sec) to aid in numerical stablitiy by simulating
% the finite time needed to vaporize
% following variable is the lag const
tc = dt/0.15; % unitless
% lagged mass vaporized (1st order lag)
MvLag = tc*(Mv-MvLag)+MvLag; % untis = kg
% sets old mass vaporized
OMv = MvLag; % units = kg
% updates old liquid nitrous mass
Mold = Mnew; % units = kg

end

