%% Basic RBC model with partial depreciation and endogenous labor
% Solved by Projection (Chebyshev Polynomials)
% Original by Jesus Fernandez-Villaverde at Haverford, July 31, 2013
% Modifiied by Kory Kantenga at University of Pennsylvania, Oct 28, 2013

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Steady State

global aalpha bbeta ddelta ppsi

aalpha = 1/3;     % Elasticity of output w.r.t. capital
bbeta  = 0.95;    % Discount factor
ddelta = 0.09;    % Rate of depreciation

% Productivity values
vProductivity = [exp(-0.04); exp(0); exp(0.04)]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000;
                 0.0041, 0.9806, 0.0153;
                 0.0000, 0.0082, 0.9837];
             
laborSteadyState = 1/3;
capitalSteadyState = laborSteadyState*((((1/bbeta)+ddelta-1)/aalpha)^...
    (1/(aalpha-1)));
outputSteadyState = (capitalSteadyState^aalpha)*(laborSteadyState^...
    (1-aalpha));
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState;
fprintf('\n');
fprintf('Steady State Values\n');
fprintf(' Output = %2.6f, Labor = %2.6f\n',outputSteadyState, laborSteadyState); 
fprintf(' Capital = %2.6f, Consumption = %2.6f\n',  capitalSteadyState,...
    consumptionSteadyState);
fprintf('\n');

%% 2. Calibration

ppsi = (1/consumptionSteadyState)*3*(1-aalpha)*((((1/bbeta)+ ddelta - 1)/...
    aalpha)^(aalpha/(aalpha-1)));
cover_k = 0.5;
klow = (1-cover_k)*capitalSteadyState;
khigh = (1+cover_k)*capitalSteadyState;
fprintf('Calibration\n');
fprintf(' Psi = %2.6f\n', ppsi); 
fprintf('\n');
fprintf('K Lower Limit\n');
fprintf(' K_underscore = %2.2f\n', klow); 
fprintf('K Upper Limit\n');
fprintf(' K_bar = %2.2f\n', khigh); 
fprintf('\n');

%% 3. Solving
options = optimset('Display','off');

vTime = zeros(2,1);

%Direct
nTime = tic;
theta0 = [0.2 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
k = col_points(klow,khigh,8);
resid = @(x) nonlinear2(x,k,vProductivity(1),vProductivity,...
    mTransition,klow,khigh);
[thetaD,fvalD] = fsolve(resid,theta0,options);
fprintf('Directly estimating with a guess\n');
display(fvalD)
display(thetaD)
vTime(1,1)=toc(nTime);

fprintf('\n');

%Iterative Method
nTime = tic;
thetaI = 0.2; %initialisation
for iOrder = 2:9
k = col_points(klow,khigh,iOrder);
resid = @(x) nonlinear2(x,k,vProductivity(2),vProductivity,...
    mTransition,klow,khigh);
[thetaI,fval] = fsolve(resid,[thetaI 0.0],options);
end
fprintf('Solving Iteratively...\n');
display(fval)
display(thetaI)
vTime(2,1)=toc(nTime);

fprintf('\n');
fprintf('The time comparison is...');
display(vTime)
fprintf('\n');

%% 4. Plots

figure;
thetaD = real(thetaD);
yD = cheby_approx(klow:0.1:khigh,thetaD,klow,khigh);
kgrid = klow:0.1:khigh;
plot(kgrid,yD,'--')
hold on
thetaI = real(thetaI);
yI = cheby_approx(klow:0.1:khigh,thetaI,klow,khigh);
kgrid = klow:0.1:khigh;
plot(kgrid,yI,'-r')
title('Labor Function for z=1','FontSize', 16)
xlabel('Capital Stock','FontSize', 12)
ylabel('Labor (Hours)','FontSize',12)
print -depsc2 HW3_LaborFunc_z1.eps

toc