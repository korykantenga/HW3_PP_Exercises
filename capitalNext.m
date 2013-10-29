function g = capitalNext(k,z,l)

global aalpha
global ddelta

%Capital next period given k, z, and l

g = -consumption(k,z,l) + z*(k^aalpha)*(l^(1-aalpha))+(1-ddelta)*k;