% HRM Sizing Script

% Johnny Jaffee

% 10/29/2018

% Rocketry at VT 2018-2019

% The following script runs a variety of simulations to help determine
% sizing for a hybrid rocket motor

clc
clear
close all

% There are 3 distinct combustion chamber sizing standards that we can size
% to, all based on AeroTech COTS 98mm motors, starting at 7680 model and
% moving up in size. For this reason, this program will consist of 3
% segments, with each exploring combinations for a specific Aerotech
% standard


% Each section will run the HRMsimFunctionWithCstarandNozEff function to
% explore how changing various parameters effect motor performance, the
% input output for the function is below:

% Nitrous blowdown model adapted from Aspire Space Technical Paper

% ---------------INPUTS-------------------------------------------------
% T = intial nitrous oxide temperature (C)
% fuelType = type of fuel, 1 = parrafin, 2 = HTPB
% chambSt = standard the combustion chamber is modelled after from aerotech
% casing sizes. This includes 1 = RMS 98/7680, 2 = RMS 98/10240, 3 = RMS 98/15360
% OFint = intial O/F (in terms of mass, not mass flow rate)
% OxTankID = Oxidizer tank ID (mm)
% cstarEff = cstar effciency
% nozEff = nozzle effcieny
% grainID = fuel grain inner diameter (mm)
% ullageSpace = %  of empty space in the oxy tank for ullage
% K = loss coeffcient for oxidizer feed system
% N = # of injectors
% injD = diameter of each injector (mm)
% dt = time step interval (sec)
% -------------------------------------------------------------------------


% --------------------OUTPUTS---------------------------------------------
% fuelOD = OD of fuel grain (mm)
% fuelLen = fuel grain length(mm)
% fuelMass = mass of fuel (g)
% nitrousMassInt = intial mass of nitrous (kg)
% TankV = nitrous tank volume (m^3)
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
% At = area of nozzle throat (m^2) (no c* correction)
% EndD = end diameter of nozzle throat (m)
% NozLen = length of nozzle (m)
% ThroatD = diamter of throat (m)
% pFuelUsed = percentage of fuel burned
% pNitrousUsed = percent nitrous used
% ---------------------------------------------------------------------


%%  --------------INTIALIZATION OF CONSTANT VARIABLES----------------------
% these variables will not be modified, and rather set constant for all RMS
% sizes
T = 24; % temp nitrous to start, C
fuelType = 1; % assuming parrafin fuel
cstarEff = 0.6; % assuming cstarEff of 0.5;
nozEff = 0.75; % assuming nozzle effciency of 0.75;
grainID = 1.125*25.4; % keeping the grain ID the same as for the solid equiv
ullageSpace = 0.2; % 25%  free space in nitrous tank
K = 2; % loss coeffcient of 2 for the feed system
N = 1; % single injector
dt = 0.01; % 0.01 sec time steps
minBurnTime = 5; % sec
dInjD = (1/16)*25.4; % value to change injD by (mm)
TankV = 15000*(1*(10^-6)); % m^3

OxTankID = 5*25.4; % mm
injDInt = (1/8)*25.4; % staring injector diameter (mm)
maxMotorLen = 6; % max allowable length of motor (ft),
minImpulseLiq = 9000; % min allowable impulse, Ns
% fprintf('The max motor length is % .1f ft\n', maxMotorLen);
% fprintf('The min motor impulse is % .1f Ns\n', minImpulseLiq);
fprintf('The intial nitrous Temp used for all calculations is % .1f C\n', T);
fprintf('The loss coeffcient K used for all calculations is % .1f\n', K);
fprintf('The intial nitrous Temp used for all calculations is % .1f C\n', T);
fprintf('The intial ullage space used for all calculations is % .1f decimal percent\n', ullageSpace);
fprintf('The intial grain ID used for all calculations is % .1f mm\n', grainID);
if fuelType == 1
    disp('parrafin fuel selected for all motor combos');
else
    disp('HTPB selected for all motor combos');
end
disp('Only thrust curves will be shown for motors that can be built within max length and meet min Impulse');
% Inputs:
% T = Temperature to be evaluated (C)

% Outputs:
% Hv = Latent heat of vaporization at given temp (J/kg)
% Cp = Specific heat capacity w/ const. pressure at given temp (j/kg)
% LiqD = liquid nitrous density at given temp (kg/m^3)
% VapD = vapor density at given T (kg/m^3)
% VapPres = Vapor Pressure at given Temp. (bar)
[ Hv,Cp,LiqD,VapD,VapPres ] = NitrousLinRegressforPhysProperties( T );


%%  ------------------END CONSTANT VARIABLE SECTION--------------------------

% intial values for variables
burnTimeInt = 20; % random high value for burn time


% Testing each RMS size will work as follows, for a given injector size, OF
% ratios will be run through. Impulse will be found for each OF, as will
% tank length be plotted as a function of tank ID for a given OF.


%%  --------------Strings For Figures--------------------------------------

strChamb = 'RMS 98/';
strCstar = num2str(cstarEff);
strNozEff = num2str(nozEff);
assumLine = strcat('cstar Effciency: ', strCstar, '| Nozzle Effciency: ', strNozEff);


%%  -----------------Section 1: RMS 98/7680----------------------------------
chambSt = 3; % RMS 98/15360 chamber
injD = injDInt; % mm
burnTime = burnTimeInt; % sec
% variable to control inner for loop from running
runInner = 1;
strMotorNum = '15360';

% keep increasing injector size until min burn time is reached
while burnTime > minBurnTime % loop 1: change injector size (increase)
    
    strInjD = num2str(injD);
    [ fuelOD, fuelLen, fuelMass, nitrousMassInt, TankLen, Time, fuelOut, Tc, impulse, Tp, VapPresp, Pcp, mdotp, Fp, nitrousMdotp, mfuelrp, totalMotorLen, impulseLiq, At, EndD, NozLen, ThroatD, pFuelUsed, pNitrousUsed  ] = HRMsimFucntionWithCstarandNozEffGivenTank(T, fuelType, chambSt, cstarEff, TankV, OxTankID, nozEff, grainID, ullageSpace, K, N, injD, dt);
    strMotor = strcat('Motor Config: ',strChamb,strMotorNum,'| Inj Diam (mm): ',strInjD);
    
    figure
    tL1 = 'Thrust vs Time';
    liqImpStr = num2str(impulseLiq);
    impLine = strcat('Total Impulse (Liquid Portion): ',liqImpStr, ' Ns');
    fUsedLine = strcat('Percentage of Fuel Burned: ', num2str(pFuelUsed));
    nUsedLine = strcat('Percentage Nitrous Used: ', num2str(pNitrousUsed));
    
    plot(Time, Fp);
    title({tL1; strMotor; assumLine; impLine; fUsedLine; nUsedLine});
    xlabel('Time (sec)');
    ylabel('Thrust (lbf)');
    grid on;
    
    
    
    injD = injD + dInjD; % mm
    burnTime = Time(end);
    
end

display(totalMotorLen);


