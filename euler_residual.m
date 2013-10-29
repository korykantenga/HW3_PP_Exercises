function y = euler_residual(k,z,vProductivity,mTransition,theta,a,b)

global aalpha ddelta bbeta

%Integrand for Residual Function (1xJ vector)

l = laborSMZ(k,theta,a,b);
Ez = zeros(1,length(vProductivity));
for j=1:length(vProductivity)
    newl = laborSMZ(capitalNext(k,z,l),theta,a,b);
    Ez(1,j) = (consumption(k,z,l)/...
        consumption(capitalNext(k,z,newl),vProductivity(j),newl))*...
        (aalpha*vProductivity(j)*...
    ((capitalNext(k,z,newl)/newl)^(aalpha-1))...
    +1-ddelta);
end

%Compute integral by Tauchen
E = Ez*mTransition((vProductivity==z),:)';
y = bbeta*E-1;