function [Tflame,qw]= flametempcombustioneff(Ti,Se,Cper,Fuel)
% James Riet, 10/16/17
% Summary
% Function takes selectivity, initial temperature, percent conversion of
% fuel and fuel type, paraffin or HTPB and outputs a flame tempera ture

% INPUTS
% Ti = Initial temperature in K
% Se = Selectivity in moles CO2 prod./moles CO prod.
% Cper = Percent conversion of fuel in percentage
% Fuel = Fuel type, either 'P' for paraffin wax or 'HTPB' for hydroxyl
% terminated polybutadiene rubber

% OUTPUTS
% Tflame = Flame temperature in degrees K

% NOTE: Assumes no heat transfer and constant pressure adibiotic flame temperature
% flame temperature calculated via mole percent CO2 and CO reaction
Hwa=23478; % enthlapy of combustion of parrafin wax with N2O
Hwb=13230; % enthalpy of reaction of N2O and wax to form CO, water, and nitrogen gas
Hpa=3109; % enthalpy of combustion of HTPB rubber with N2O
Hpb=1645; % enthalpy of reaction of N2O and HTPB to form CO, water, and nitrogen gas
CpCOO=.05623; % Constant pressure heat capacities at 1200K of every product that will be present in at least one of the reactions,
CpCO=.03416;  % 1200K used as intermediate heat capacity to model the average heat capacity of the products over a wide range of temperatures
CpHHO=.04365;
CpNN=.0337;
CpNNO=.03872;
if Fuel == 'P'
    wnCOO=(Se/(Se+1))*28*(Cper/100); % moles of each product formed given input parameters per one mole of fuel
    wnCO=(1/(Se+1))*28*(Cper/100);
    wnHHO=(Cper/100)*29;
    wnNN=(Cper/100)*(Se/(1+Se))*85 + (Cper/100)*(1/(Se+1))*57;
    wnNNO=(1-(Cper/100))*85;
    qw=((Hwa*(Se/(Se+1))*(Cper/100)) + (Hwb*(1/(Se+1)*(Cper/100))));
    DeltaT =  qw/(wnCOO*CpCOO + wnCO*CpCO+ wnHHO*CpHHO + wnNN*CpNN + wnNNO*CpNNO); % energy balance that calculates the delta H of reaction then balances that with products formed
    Tflame=Ti+DeltaT; % find flame temp from inital temp
end
if Fuel == 'HTPB'
    pnCOO=(Se/(Se+1))*4*(Cper/100); % moles +of each product formed given input parameters per one mole of fuel
    pnCO=(1/(Se+1))*4*(Cper/100);
    pnHHO=(Cper/100)*3;
    pnNN=(Cper/100)*(Se/(1+Se))*11 + (Cper/100)*(1/(Se+1))*7;
    pnNNO=(1-(Cper/100))*85;
    qp=((Hpa*(Se/(Se+1))*(Cper/100)) + (Hpb*(1/(Se+1))*(Cper/100)));
    DeltaT = qp/(pnCOO*CpCOO + pnCO*CpCO+ pnHHO*CpHHO + pnNN*CpNN + pnNNO*CpNNO);
    Tflame=Ti+DeltaT; % find flame temp from inital temp
end

% fprintf('\nThe flame temperature is % 7.1f K\n',Tflame)
end