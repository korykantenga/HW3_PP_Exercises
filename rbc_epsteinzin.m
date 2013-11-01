// Basic RBC Model with CRRA preferences
//
// See basic_notes.tex for details
//
// Jesus Fernandez-Villaverde
// Philadelphia, April 9, 2009
// Modified by Kory Kantenga (University of Pennsylvania, Oct 2013)

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables (20=1+6+2+2+4+5)
//----------------------------------------------------------------
 
var 

// Utility (1)
v

// Allocation (6)
y c i l k g

// Prices (2)
r w

// Shadow Prices (2)
lambda m

// Productivity, Taxes & Preference Shock (4)
z taul tauk eta

//Auxillary Variables (5)
AuxI AuxK X Y U;


//----------------------------------------------------------------
// 2. Exogenous variables (4)
//----------------------------------------------------------------
 
varexo 
ez;

varexo
etaul;

varexo

etauk;

varexo
eeta;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Preferences
bbeta zzeta kkappa 

// Technology
aalpha ddelta rrhoz rrhoeta

//Taxes
rrhotau1 rrhotauk

// S.D.'s stochastic processes
ssigmaz
ssigmazlag
ssigmaeta
sigmataul
sigmatauk;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences
bbeta      = 0.99;
zzeta      = 0.5;
kkappa     = -9;

// Technology
aalpha     = 0.33;
ddelta     = 0.03;
rrhoz      = 0.95; 
rrhoeta    = 0.95;
rrhotau1   = 0.95; 
rrhotauk   = 0.90;

// Stochastic processes
ssigmaz    = 0.007;
ssigmazl   = 0.001;
ssigmaeta  = 0.002;
sigmataul  = 0.005;
sigmatauk  = 0.008;

//----------------------------------------------------------------
// 5. Computation of Steady State
//----------------------------------------------------------------

// AR Processes
eta_ss  = 1;
taul_ss = 0.25;
tauk_ss = 0.35;

// Endogenous Variables
r_ss  = (1/bbeta-1+ddelta)/(1-tauk_ss);
w_ss  =(1-aalpha)*(aalpha/r_ss)^(aalpha/(1-aalpha));
l_ss  = sqrt((1-taul_ss)*w_ss/(((aalpha/r_ss)^(aalpha/(1-aalpha))*(1-(1-aalpha)*taul_ss-aalpha*tauk_ss)-ddelta*(aalpha/r_ss)^(1/(1-aalpha)))));
k_ss  = ((r_ss/(aalpha*(l_ss^(1-aalpha))))^(1/(aalpha-1)));
i_ss  = ddelta*k_ss;
y_ss  = (k_ss^aalpha)*(l_ss^(1-aalpha));
w_ss  = (1-aalpha)*y_ss/l_ss;
c_ss  = (1-taul_ss)*w_ss*l_ss + (1-tauk_ss)*r_ss*k_ss - i_ss;
g_ss  = taul_ss*w_ss*l_ss + tauk_ss*r_ss*k_ss;
v_ss = (log(c_ss) - ((l_ss^2)/2))/((1-bbeta)^2);
U_ss = v_ss^kkappa;
lambda_ss = real((v_ss^zzeta)*(1/c_ss)*((log(c_ss)- ((l_ss^2)/2))^(zzeta-1)));
m_ss = lambda_ss;
Y_ss = lambda_ss*(1-tauk_ss)*r_ss + (1-ddelta)*m_ss;
AuxK_ss = (v_ss^-10)*Y_ss;

//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model; 

  // 1. Productivity process 
  z = rrhoz*z(-1) + (ssigmaz*ez) + (0.001*ez(-1));

  // 2. Labor Tax process
  taul = 0.25 + 0.95*(taul(-1)-0.25)+0.005*etaul;

  // 3. Capital Tax process
  tauk = 0.35 + 0.9*(tauk(-1)-0.35)+0.008*etauk;

  // 4. Preference Shock process
  eta = 1+ 0.95*(eta(-1)-0.35)+0.002*eeta;

  // 5. Resource Constraint
  c = exp(z)*(k(-1)^aalpha)*(l^(1-aalpha))-(i + g);

  // 6. Goverment Budget Constraint
  g = tauk*w*l + taul*r*k(-1);

  // 7. Production function
  y = (exp(z)*(k(-1)^aalpha)*(l)^(1-aalpha));

  // 8. Wage Equation
  w = (1-aalpha)*y/l;

  // 9. Rental Rate Equation
  r = aalpha*y/k(-1);

 // 10. Capital Law of Motion
 k = (1-ddelta)*k(-1) + i*(1-(0.05*(((i/i(-1))-1)^2)));

 // 11. Intertemporal Optimality (Labor Supply Decision)
 l = ((1-taul)*w)/(c*eta);

 // 12. Auxillary Value Function
 U = v(+1)^kkappa;

 // 13. Value Function
 v = (((1-bbeta)*(log(c)-eta*((l^2)/2))^zzeta+bbeta*U^(zzeta/kkappa))^(1/zzeta));

 // 14. FOC(c)
lambda = (v^zzeta)*(1/c)*((log(c)-eta*((l^2)/2))^(zzeta-1));

 // 15. EC(i)
X = -0.1*m(+1)*((i(+1)/i)-1)*((i(+1)/i)^2);

// 16. EC(k)
Y = lambda(+1)*(1-tauk(+1))*r(+1) + (1-ddelta)*m(+1);

// 17. Auxillary I
AuxI = (v(+1)^-10)*X;

// 18. Auxillary K
AuxK = (v(+1)^-10)*Y;

// 19. FOC(i)
AuxI = ((bbeta*(U^((zzeta/kkappa)-1))*(v^zzeta))^-1)*(lambda + m*(1-0.05*(((i/i(-1))-1)*(((i/i(-1))-1)+((2*i)/i(-1))))));

// 20. Derivative of Lagrangian w.r.t. k
m = bbeta*(U^((zzeta/kkappa)-1))*AuxK*(v^zzeta);

end;

initval;
  v = v_ss;
  y = y_ss;
  c = c_ss;
  i = i_ss;
  l = l_ss;
  k = k_ss;
  r = r_ss;
  w = w_ss;
  U = U_ss;
  g = g_ss;
  lambda = lambda_ss;
  m = m_ss;
  AuxI = 0;
  X = 0;
  AuxK = AuxK_ss;
  Y = Y_ss;
  z = 0;
  taul = 0.25;
  tauk = 0.35;
  eta = 1;
  ez = 0;
  etaul = 0;
  etauk = 0;
  eeta  = 0;
end;

shocks;
  var ez    = 1;
  var etaul = 1;
  var etauk = 1;
  var eeta  = 1;
end;

steady;

stoch_siml(hp_filter = 1600, irf = 20, order = 3);