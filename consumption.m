function c = consumption(k,z,l)

%Determine consumption given k, z, and l(k,z)
global aalpha
global ppsi

c1 = ((1-aalpha)*z*(k/l)^(aalpha))/(ppsi*l);
c=max(c1,0.0000001);