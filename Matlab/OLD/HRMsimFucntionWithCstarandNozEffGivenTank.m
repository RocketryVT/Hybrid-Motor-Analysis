function [ fuelOD, fuelLen, fuelMass, nitrousMassInt, TankLen, Time, fuelOut, Tc, impulse, Tp, VapPresp, Pcp, mdotp, Fp, nitrousMdotp, mfuelrp, totalMotorLen, impulseLiq, At, EndD, NozLen, ThroatD, pFuelUsed, pNitrousUsed  ] = HRMsimFucntionWithCstarandNozEffGivenTank(T, fuelType, chambSt, cstarEff, TankV, OxTankID, nozEff, grainID, ullageSpace, K, N, injD, dt)
% This function models a hybrid rocket motor firing, while accounting for
% c* effciency, applied to combustion chambe rpressure and Temp
%   Assumes using nitrous oxide
%   Fuel is either parrafin or HTPB,
%


% Nitrous blowdown model adapted from Aspire Space Technical Paper

%% ---------------INPUTS-------------------------------------------------
% T = intial nitrous oxide temperature (C)
% fuelType = type of fuel, 1 = parrafin, 2 = HTPB
% chambSt = standard the combustion chamber is modelled after from aerotech
% casing sizes. This includes 1 = RMS 98/7680, 2 = RMS 98/10240, 3 = RMS 98/15360
% TankV = nitrous tank volume (m^3)
% OxTankID = Oxidizer tank ID (mm)
% cstarEff = cstar effciency
% nozEff = nozzle effcieny
% grainID = fuel grain inner diameter (mm)
% ullageSpace = decimal % of empty space in the oxy tank for ullage
% K = loss coeffcient for oxidizer feed system
% N = # of injectors
% injD = diameter of each injector (mm)
% dt = time step interval (sec)

%% -------------------------------------------------------------------------


%% --------------------OUTPUTS---------------------------------------------
% fuelOD = OD of fuel grain (mm)
% fuelLen = fuel grain length(mm)
% fuelMass = mass of fuel (g)
% nitrousMassInt = intial mass of nitrous (kg)

% TankLen = length of nitrous oxide tank (m)
% Time = vector of time after ignition (sec)
% fuelOut = boolean value for if fuel ran out (0 = no, 1 = yes)
% Tc = combustion chamber flame temp (K) with cstar applied
% impulse = total motor impulse (Ns)
% Tp = vector of nitrous oxide temp (C)
% VapPresp = vector of nitrous oxide vapor pressure(bar)
% Pcp = vector of combustion chamber pressure (bar)
% mdotp = total propellant mass flow rate vector (kg/s)
% Fp = thrust vector (lbf)
% nitrousMdotp = vector of nitrous oxide mass flow rate (kg/s)
% mfuelrp = vector of fuel mass flow rate (kg/s)
% totalMotorLen = estimate of total motor length (ft)
% impulseLiq = impulse for liquid portion Ns
% At = area of nozzle throat (m^2) no c*correction
% EndD = end diameter of nozzle throat (m)
% NozLen = length of nozzle (m)
% ThroatD = diamter of throat (m)
% pFuelUsed = percentage of fuel burned
% pNitrousUsed = percent nitrous used
%% -------------------------------------------------------------------------



% Intialization-------------------------------------------------------
% begin with intial guess of vaporized nitrous mass
Mv = .00000001; % Units = Kg
% variable for atomsphereic pressure
Pa = 1.013529; % units = bar



% Determine some characteristics regarding the fuel
if fuelType == 1
    disp('Wax Fuel Selected');
    fuelDen = 0.9; % g/cm^3
    a = 0.8;
    b = 0.5036;
    Mav = 29; % g/mol
    Ftype = 'P';
end
if fuelType == 2
    disp('HTPB Fuel Selected');
    fuelDen = 0.93; % g/cm^3
    a = 0.15;
    b = 0.5;
    Ftype = 'HTPB';
    Mav = 29.9; % g/mol
end

% find nitrous properties
% value for Runi in J/kmol/k
Runi = 8314.4621; % UNITS = J/kmol/k
% ambient atmospheric pressure (bar)
pa = 1.0135; % units = bar. Pressure SL
% gets nitrous physical properties
% Inputs:
% T = Temperature to be evaluated (C)

% Outputs:
% Hv = Latent heat of vaporization at given temp (J/kg)
% Cp = Specific heat capacity w/ const. pressure at given temp (j/kg)
% LiqD = liquid nitrous density at given temp (kg/m^3)
% VapD = vapor density at given T (kg/m^3)
% VapPres = Vapor Pressure at given Temp. (bar)
[ Hv,Cp,LiqD,VapD,VapPres ] = NitrousLinRegressforPhysProperties( T );
% fuel grain OD
fuelOD = 3.375*25.4; % mm, sized for RMS 98mm size
% Find mass of fuel, nitrous, and tank size
if chambSt == 1
    disp('RMS 98/7680 Chamber');
    % fuel grain length
    fuelLen = 18.75*25.4; % mm
end
if chambSt == 2
    disp('RMS 98/10240 Chamber');
    % fuel grain length
    fuelLen = 24.813*25.4; % mm
end
if chambSt == 3
    disp('RMS 98/15360 Chamber');
    % fuel grain length
    fuelLen = 36.938*25.4; % mm
end



% volume of fuel, note that the OD and ID are being converted to radius
% and cm at same time
fuelVol = (pi*((fuelOD/20)^2)*(fuelLen/10))-(pi*((grainID/20)^2)*(fuelLen/10)); % cm^3
fuelMass = fuelDen*fuelVol; % grams
nitrousMassInt = (TankV-(TankV*ullageSpace))*LiqD;
% inital mass of liquid nitrous
Mold = nitrousMassInt; % kg

% calculate the volume of nitrous oxide
nitrousVol = Mold/LiqD; % m^3
% Find Tank volume, with ullage applied

% calculate the length of the tank
TankLen = TankV/(pi*((OxTankID/2000)^2)); % m
% intial chamber pressure, assume half of ox tank
Pc = VapPres*0.5; % bar
% specific heat ratio
k=1.4; % UNITLESS
% universal gas constant in kj/mol
runi = 0.0083145; % units = kj/mol
% variable for loss coeffcient. Empirical. Default =2
K = 2; % unitless
% calculates the area of each injector
injA = (((injD/1000)*0.5)^2)*pi; % units = m^2
% expression for the denominator on the mass flow rate eq from aspire space
Dloss = (K/((N*injA)^2))*(1/100000); % units = 1/m^4
% calculates the change in pressure from manifold to post injector
dPres = VapPres-Pc; % units = bar
% converts the pressure difference to units of bar
% dPres = (dPres*0.01); % units = bar
% calculates mass flow rate of nitrous
nitrousMdot = sqrt((2*LiqD*dPres)/Dloss); % units = kg/sec
% fuel grain length and diameter in cm
Cl = fuelLen/10; % cm
Cd = fuelOD/10; % cm
% fuel grain port radius in cm
portr = grainID/20; % cm
% thickness of fuel in cm
fuelThick = (Cd/2)-portr; % cm
% Outputs for function called below
% mfuelr = fuel mass flow rat % units =kg/sec
% r = regression rate % units = mm/sec
% mdot = total massflow kg/s
% starts with dt = 0 as to not include erosion effects


[ mfuelr,r,mdot,portr ] = FuelFlow_RegressRate( portr,Cd, Cl, nitrousMdot, a, b,fuelDen,0);
% INPUTS
% Ti = Initial temperature in K
% Se = Selectivity in moles CO2 prod./moles CO prod.
% Cper = Percent conversion of fuel in percentage
% Fuel = Fuel type, either 'P' for paraffin wax or 'HTPB' for hydroxyl
% terminated polybutadiene rubber

% OUTPUTS
% Tc = Flame temperature in degrees K
% accounting for cstarEff
[Tc,qw]= flametempcombustioneffFinal((70+273),(15),(95),Ftype);
% finds the temp (k) at the throat of the nozzel
Tt = (Tc/(1+((k-1)/2))); % UNITS = KELVIN
% At = nozzle throat area % units = m^2, no c* correction
[ At,EndD,NozLen,ThroatD ] = NozDimen( mdot,Pc,Mav,Tt,Pa ); % units = m^2
% variable for totaltime
Time = 0; % units = sec
TimeLast = -dt;
% total nitrous mass
Mtot = Mold+Mv; % units = kg

% find the intial Mnew
Mnew = (TankV-((Mtot/VapD)))/((1/LiqD)-(1/VapD)); % units = kg
% calculates the mass of nitrous vapor
MassVap = Mtot-Mnew; % units = kg
% resets the intial Mold after figuring out how much of the nitrous is
% vapour
Mold = Mtot-MassVap; % units = kg
Mnew = Mold-0.0001;% units =kg
% sets OMv = Mv
OMv = Mv; % units = kg
% sets MvLag = OMv;
MvLag = OMv; % units = kg
% sets Mnew = Mold;
% Mnew = Mold-.000001; % units = kg
% loop counter
lpCount = 1;
% nitrousMdot = 0;
PcLag=Pc;
MdotLag = mdot;
% ^Intial conditions----------------------------------------------------

% Variable for impluse (Ns)
impulse = 0; % Ns

% boolean variable to determine what runs out first, oxidizer or fuel, 1 =
% fuel ran out first, 0 = oxidizer
fuelOut = 0;

% INPUTS
% Ti = Initial temperature in K
% Se = Selectivity in moles CO2 prod./moles CO prod.
% Cper = Percent conversion of fuel in percentage
% Fuel = Fuel type, either 'P' for paraffin wax or 'HTPB' for hydroxyl
% terminated polybutadiene rubber

% OUTPUTS
% Tc = Flame temperature in degrees K
% accounting for cstarEff
[Tc,qw]= flametempcombustioneffFinal((70+273),(15),(95*cstarEff),Ftype);
% finds the temp (k) at the throat of the nozzel
Tt = (Tc/(1+((k-1)/2))); % UNITS = KELVIN



% Models Liquid emptying out of nitrous tank, does not include vapor phase
% Mnew=0;
while ((Mnew > 0.15) && (fuelOut == 0))
    % Mnew > 1 && Mold > 1
    [ Pc ] = ChamberPressure( mdot,At,Tt,Mav, k ); % bar
    
    % apply cstarEff correction
    Pc = Pc*cstarEff; % bar
    tc = dt/0.15; % unitless
    % lagged pc (1st order lag)
    PcLag = tc*(Pc-PcLag)+PcLag; % untis = bar
    Pc=PcLag;
    % Pc = 0.7*VapPres;
    % the function below simulates the liquid nitrous emptying out in given
    % time steps
    [ T,Mv,Mold,Mtot,Mnew,OMv,MvLag,nitrousMdot,MassVap,VapPres ] = LiqNitrousEmpty( T,Mv,Mold,TankV,Pc,N,injD,dt,Mtot,Mnew,OMv,MvLag,nitrousMdot,MassVap );
    nitrousRemain = Mnew; % kg
    % calculates the pressure drop from injector to chamber
    Pcdrop = ((VapPres-Pc)/VapPres)*100; % units = %
    % checks to see if within safe level (20% )
    if Pcdrop < 20
        disp('Warning: Unsafe Pressure Drop Across Injectors');
        display(Pcdrop);
    end
    % disp((Cd/2)-portr);
    [ mfuelr,r,mdot,portr] = FuelFlow_RegressRate( portr,Cd, Cl, nitrousMdot, a, b,fuelDen,dt);
    % check to see is the fuel is all regressed
    if (portr >= ((Cd/2)-((Cd/2)*0.05)))
        fuelOut = 1;
        disp('Out of Fuel');
    end
    % lagged mdot (1st order lag)
    MdotLag = tc*(mdot-MdotLag)+MdotLag; % untis = kg/s
    mdot=MdotLag;
    % exhaust velocity calculation NOT ADJUSTED FOR EQUIVALENT PRESSURE. PA
    % SHOULD = PE
    Veq = (sqrt(((2*k)/(k-1))*((Runi*Tc)/Mav)*(1-((pa/Pc)^((k-1)/k))))); % UNITS = m/s
    % apply nozzle efffciency correction
    Veq = Veq*nozEff;
    % thrust
    F = (Veq*mdot)/4.45; % UNITS = Pounds
    impulse = impulse + (F*4.45*dt); % Ns
    % [ Pc ] = ChamberPressure( mdot,At,Tt,Mav );
    if Mnew > Mold
        Mnew = 0;
    end
    
    
    
    
    
    Mvp(lpCount) = Mv;
    % MassLostp(lpCount+1) = MassLost;
    Moldp(lpCount) = Mold;
    Mnewp (lpCount) = Mnew;
    Mtotp(lpCount) = Mtot;
    Tp(lpCount) = T;
    % dTp(lpCount+1) = dT;
    VapPresp(lpCount)=VapPres;
    Pcp(lpCount) = Pc;
    mdotp(lpCount) = mdot;
    Fp(lpCount) = F;
    nitrousMdotp(lpCount) = nitrousMdot;
    mfuelrp(lpCount) = mfuelr; % kg/s
    % updeates time elapsed
    
    Time(lpCount) = TimeLast + dt; % units = sec;
    TimeLast = Time(lpCount); % sec
    lpCount = lpCount+1;
end
fprintf('Ending (Liq. Phase) Pressure Drop Across Injector = % .3f Percent\n',Pcdrop);
impulseLiq = impulse; % Ns

% initalizes variables to be used in simulating the vapor portion of the
% tank emptying
MassVap = Mtot; % units = kg
% intializes intial compressibility factor
[ Zint ] = NitrousCompressFact( VapPres ); % unitless
% mass of vapor
MassVapint = MassVap; % units = kg;
% intializes intial Temp
Tint = T; % units = kg
% resets nitrousMdot
nitrousMdot = 0; % units = kg/sec
% store intial vapor pressure
VapPresint = VapPres; % units = bar
% store intial vapor density
VapDint = VapD; % units = kg/m^3
Pc = VapPres*0.6;
c = 1;
Tme = 0;
Tmes = 0;
PcLag = Pc;
VapLag = VapPres;



while ((MassVap > 0.1) && (VapPres > 4.5) && (fuelOut == 0))
    tc = dt/0.15; % unitless
    [ Pc ] = ChamberPressure( mdot,At,Tt,Mav, k );
    % apply cstar effciency correction
    Pc = Pc*cstarEff; % bar
    % lagged pc (1st order lag)
    PcLag = tc*(Pc-PcLag)+PcLag; % untis = bar
    Pc=PcLag;
    [ VapD,nitrousMdot,T,MassVap,VapPres] = NitrousGasEmptyFunct( Pc,Zint,MassVapint,Dloss,VapD,nitrousMdot,dt,Tint,T,MassVap,VapPres,VapPresint,VapDint );
    nitrousRemain = MassVap; % kg
    % Pc = VapPres*0.6;
    VapLag = tc*(VapPres-VapLag)+VapLag; % untis = bar
    VapPres=VapLag;
    [ mfuelr,r,mdot,portr] = FuelFlow_RegressRate( portr,Cd, Cl, nitrousMdot, a, b,fuelDen,dt);
    % checks that there is still fuel left
    if (portr*2) >= (Cd-(Cd*0.05))
        %  MassVap = 0; % units = kg
        disp('Out of Fuel Before Oxidizer');
        fuelOut = 1;
        
    end
    % disp((Cd/2)-portr);
    % lagged mdot (1st order lag)
    MdotLag = tc*(mdot-MdotLag)+MdotLag; % untis = kg/s
    mdot=MdotLag;
    % exhaust velocity calculation NOT ADJUSTED FOR EQUIVALENT PRESSURE. PA
    % SHOULD = PE
    Veq = (sqrt(((2*k)/(k-1))*((Runi*Tc)/Mav)*(1-((pa/Pc)^((k-1)/k))))); % UNITS = m/s
    % apply nozzle effciecy correction
    Veq = Veq*nozEff;
    % thrust
    F = (Veq*mdot)/4.45; % UNITS = Pounds
    % calc impulse
    impulse = impulse+(F*4.45*dt); % Ns
    
    VapDp(c) = VapD;
    VapPrespp(c) = VapPres;
    Tme(c) = Tmes+dt;
    Tmes = Tme(c);
    VapPresp(lpCount)=VapPres;
    Pcp(lpCount) = Pc;
    Time(lpCount) = TimeLast + dt; % units = sec;
    TimeLast = Time(lpCount);
    Fp(lpCount) = F;
    nitrousMdotp(lpCount) = nitrousMdot;
    mfuelrp(lpCount) = mfuelr; % kg/s
    c = c+1;
    lpCount = lpCount+1;
    
end

if (fuelOut == 0)
    disp('Oxidizer out before fuel');
end
fprintf('The total impulse is % .1f N*s\n', impulse);

totalMotorLen = (((fuelLen/10)/2.54)/12)+(((TankLen*100)/2.54)/12); % ft
if (fuelOut == 1)
    pFuelUsed = 100;
else
    
    fuelVolFinal = (pi*((fuelOD/20)^2)*(fuelLen/10))-(pi*((portr)^2)*(fuelLen/10)); % cm^3
    fuelMassFinal = fuelDen*fuelVolFinal; % grams
    fuelBurned = fuelMass-fuelMassFinal; % g
    pFuelUsed = (fuelBurned/fuelMass)*100; % percent
end
pNitrousUsed = ((nitrousMassInt-nitrousRemain)/nitrousMassInt)*100; % percent
end

