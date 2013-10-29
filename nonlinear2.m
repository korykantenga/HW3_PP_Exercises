function y = nonlinear2(theta,k,z,vProductivity,mTransition,a,b)

% nonlinear system of equations

y = zeros(length(k),1);

for i = 1:length(k)
    y(i,1) = ...
        euler_residual(k(i),z,vProductivity,mTransition,theta,a,b);
    
end