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
theta0 = [0.3333 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
k = col_points(klow,khigh,8);
resid = @(x) nonlinear2(x,k,vProductivity(2),vProductivity,...
    mTransition,klow,khigh);
[thetaD,fvalD] = fsolve(resid,theta0,options);
fprintf('Directly estimating with a guess\n');
display(fvalD)
display(thetaD)
vTime(1,1)=toc(nTime);

fprintf('\n');

%Iterative Method
nTime = tic;
thetaI = 0.3333; %initialisation
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

% Entire labor function
theta0 = [0.3333 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
k = col_points(klow,khigh,8);
[theta1,fval1] = fsolve(@(x) ...
    nonlinear2(x,k,vProductivity(1),vProductivity,...
    mTransition,klow,khigh),theta0,options);
[theta2,fval2] = fsolve(@(x) ...
    nonlinear2(x,k,vProductivity(2),vProductivity,...
    mTransition,klow,khigh),theta0,options);
[theta3,fval3] = fsolve(@(x) ...
    nonlinear2(x,k,vProductivity(3),vProductivity,...
    mTransition,klow,khigh),theta0,options);


%% 4. Plots

kgrid = klow:0.1:khigh;

figure;
thetaD = real(thetaD);
yD = cheby_approx(klow:0.1:khigh,thetaD,klow,khigh);

plot(kgrid,yD,'-o')
hold on
thetaI = real(thetaI);
yI = cheby_approx(klow:0.1:khigh,thetaI,klow,khigh);
plot(kgrid,yI,'-r')
title('Labor Function for z=0','FontSize', 16)
xlabel('Capital Stock','FontSize', 12)
ylabel('Labor (Hours)','FontSize',12)
xlim([klow khigh])
legend('Direct Guess','Iterative Guess','Location','Best')
print -depsc2 HW3_LaborFunc_z1.eps

figure;
theta1 = real(theta1);
y1 = cheby_approx(klow:0.1:khigh,theta1,klow,khigh);
plot(kgrid,y1)
hold on;
theta2 = real(theta2);
y2 = cheby_approx(klow:0.1:khigh,theta2,klow,khigh);
plot(kgrid,y2,'-r')
hold on;
theta3 = real(theta3);
y3 = cheby_approx(klow:0.1:khigh,theta3,klow,khigh);
plot(kgrid,y3,'-g')
hold on;
title('Labor Policy Function','FontSize', 16)
xlabel('Capital Stock','FontSize', 12)
ylabel('Labor (Hours)','FontSize',12)
xlim([klow khigh])
legend('z=-0.04','z=0','z=0.04','Location','Best')
print -depsc2 HW3_LaborFunc.eps


%% 5. Recover and Plot Consumption and Capital Policy Functions
laborPolicy = [y1;y2;y3]';

consumptionPolicy = zeros(length(kgrid),length(vProductivity));
capitalPolicy = zeros(length(kgrid),length(vProductivity));
for i=1:length(kgrid)
    for j=1:length(vProductivity)
        consumptionPolicy(i,j)=consumption(kgrid(i),...
            vProductivity(j),laborPolicy(i,j));
        capitalPolicy(i,j)=capitalNext(kgrid(i),...
            vProductivity(j),laborPolicy(i,j));
    end
end
figure;
plot(kgrid,consumptionPolicy)
title('Consumption Policy Function','FontSize', 16)
xlabel('Capital Stock','FontSize', 12)
ylabel('Consumption Units','FontSize',12)
xlim([klow khigh])
legend('z=-0.04','z=0','z=0.04','Location','Best')
print -depsc2 HW3_ConsFunc.eps
figure;
plot(kgrid,capitalPolicy)
title('Capital Policy Function','FontSize', 16)
xlabel('Capital Stock','FontSize', 12)
ylabel('Capital','FontSize',12)
xlim([klow khigh])
legend('z=-0.04','z=0','z=0.04','Location','Best')
print -depsc2 HW3_CapitalFunc.eps

toc