function [ VapD,nitrousMdot,T,MassVap,VapPres] = NitrousGasEmptyFunct( Pc,Zint,MassVapint,Dloss,VapD,nitrousMdot,dt,Tint,T,MassVap,VapPres,VapPresint,VapDint )
% the following function simulates nitrous vapor emptying form a tank
% Developed using the aid of Aspire Space technical paper on the modeling
% of nitrous tank emptying


% Inputs
% Pc = combustion chamber pressure % units = bar
% Zint = intial compressiability factor % unitless
% MassVapint = inital mass of nitrous vapor % units = kg
% Dloss = pressure loss coeffcient for mass flow rate calculations
% dt = time step interval % units = sec
% Tint = inital Temp % units = C
% T = current Temp % units = C
% MassVap = current mass of nitrous vapor % units = kg
% VapD = Vapor Density % units = kg/m^3
% Outputs
% VapD = nitrous vapor density % units = kg/m^3
% nitrousMdot = nitrous mas flow rate % units = kg/sec
% T = current temperature of nitrous % units = C
% MassVap = mass of nitrous vapor 5units = kg/sec





% sets OnitrousMdot = nitrousMdot
OnitrousMdot = nitrousMdot; % units = kg/s
% int = intial
% calculates pressure difference between chamber and Nitrous
dPres = VapPres-Pc; % units = bar
% dPres=10;
% calculates nitrous mdot
nitrousMdot = sqrt(((2*VapD*dPres)/Dloss)); % units = kg/sec
% calculates the mass lost via Addams Second Order Integration forumla:
% Xn = X(n-1)+dt/2 * ((3*Xdot(n-1)-Xdot(n-2))
MassLost = 0.5*dt*(3*nitrousMdot-OnitrousMdot); % units = kg
% calculate the new total mass in the tank
MassVap = MassVap-MassLost; % units = kg
% variable for specific heatratio of nitrous
y = 1.3; % unitless
% guesses Z by calling the function below
% Inputs
% VapPres = Pressure to find compressibility factor at. % units = bar
% Outputs
% Z = compressiblity factor at given pressure % unitless
[ Zguess ] = NitrousCompressFact( VapPres ); % unitless
% multiplier/divider factor to adjust Zguess value
% Zchange = 1.11;% unitless

% variable for error percent
ErrorPercent = 100; % units = % 
% variable for desired error marging percentage
ErrTol = 5; % units = % 
% variable for step size
step = 1/0.9; % unitless
% refines loop to get better convergence
Oaim = 2; % unitless
aim = 1; % unitless
% intializes Z such that it it much different than Zguess to get loop to
% run
Z = Zguess*2; % unitless
% loopcount variable
count = 0; % unitless
% loop to find Z to a given error tolerance
while  ErrorPercent > ErrTol
    % ErrorPercent > ErrTol && count < 1000
    % calculates Temp based on guess
    T = Tint*(((Zguess*MassVap)/(MassVapint*Zint))^(y-1)); % units = C
    % calculates pressure based on guess
    VapPres = VapPresint*((T/Tint)^(y/(y-1))); % units = bar
    % calulates new Z value
    [ Z ] = NitrousCompressFact( VapPres ); % unitless
    % resets old aim
    Oaim = aim; % unitless
    
    % Increments the Zguess
    if (Zguess < Z)
        % increments up
        Zguess = Zguess*step; % unitless
        aim =1;
    end
    if Zguess > Z
        Zguess = Zguess/step; % unitless
        aim = -1;
    end
    % checks to see if target was overshot, and if so reduces the step
    % size
    if aim == -Oaim
        step = sqrt(step); % unitless
    end
    
    % calculates error percent
    ErrorPercent = ((abs((Zguess-Z)/Z))*100); % units = % 
    % display(ErrorPercent);
    % increments count
    count = count+1;
    % increases error tolerance if program takes too long
    if count > 1000
        ErrTol = ErrTol+(ErrTol*0.1); % increases error tol by 10% 
    end
end
VapD = VapDint*((T/Tint)^(1/(y-1))); % units = kg/m^3
% calculates pressure based on guess
% VapPres = VapPresint*((T/Tint)^(y/(y-1))); % units = bar




end

