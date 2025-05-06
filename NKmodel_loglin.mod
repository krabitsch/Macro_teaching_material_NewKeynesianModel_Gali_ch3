// NKmodel_loglin.mod
// codes accompanying courses Advanced Macro II, Macro Models and Methods, MA Economics, WU Vienna, Katrin Rabitsch
//
// Solving the baseline NK model (as in Gali (2008), chapter 3) with Dynare, model variables specified in log-linear terms

var gap pi i v a y n m w y_bar r deltam y_over_n i_ann r_ann pi_ann;

varexo ev ea;

parameters betta sigma varphi mu xsi theta alpha lambda kappa phi_p phi_x rho_a rho_v sigma_a sigma_v y_st ybar_st n_st w_st pi_st gap_st i_st r_st m_st deltam_st v_st a_st;

// households 
betta  = 0.99;                                                        // this value implies a riskless annual return of about 4% 
sigma  = 1;                                                           // relative risk aversion coefficient 
varphi = 5;                                                           // inverse of the labour supply elasticity
mu     = 4;                                                           // interest semi-elasticity of money demand

// firms
xsi    = 0.75;                                                        // firms not adjusting the prices in each period 
theta  = 9;                                                           // elasticity of substitution across goods 
alpha  = 1/4;                                                         // when = 0, constant return technology
lambda = (1-betta*xsi)*(1-xsi)/xsi;
kappa  = lambda*((1-alpha)*sigma+varphi+alpha)/(1-alpha+alpha*theta);

// monetary policy
phi_p  = 1.5;
phi_x  = 0.5/4;

// exogenous processes
rho_a   = 0.9; 
rho_v   = 0.5; 
sigma_a = 0.007; 
sigma_v = 0.0025;    


model(linear);

    // Dynamic IS equation
    gap = gap(+1) - (1/sigma) * (i - pi(+1)  - sigma*(1+varphi)/((1-alpha)*sigma+varphi+alpha) * (a(+1)-a) );

    // New Keynesian Phillips Curve (NKPC)
    pi = betta * pi(+1) + kappa * gap;

    // Monetary Policy: Taylor Rule 
    i =  phi_p * pi + phi_x * gap + v;

    // Exogenous Process for Monetary Shock
    v = rho_v * v(-1) + ev; 

    // Exogenous Process for TFP Shock
    a = rho_a * a(-1) + ea;


    // Other variables carried in the system (not required to obtain solution)
    
    // labour supply 
    w = sigma * y + varphi * n;

    // money demand equation
    m = y - mu * i;
    
    //production function
    y = a + (1-alpha) * n;

    // natural (flex price) level of output
    y_bar = (1+varphi)/((1-alpha)*sigma+varphi+alpha) * a; 

    // output under sticky prices
    y = gap + y_bar;
    
    // money growth
    deltam = m - m(-1) + pi; 

    // real interest rate
    r = i - pi(+1);

    // labor productivity
    y_over_n = y - n; 

    // Annualized nominal interest rate
    i_ann=4*i;

    // Annualized real interest rate
    r_ann=4*r;

    // Annualized inflation
    pi_ann=4*pi;

end;

initval;
  gap      = 0;
  pi       = 0;
  i        = 0;
  v        = 0;
  a        = 0;
  y        = 0;
  n        = 0;
  m        = 0;
  w        = 0;
  y_bar    = 0;
  r        = 0; 
  deltam   = 0;
  y_over_n = 0;
  i_ann    = 0;
  r_ann    = 0;
  pi_ann   = 0;

end;

steady;

check;

shocks;
// for simulation with appropriate shock s.d.:
//var ea = sigma_a^2; 
//var ev = sigma_v^2;
// to obtain Dynare auto-generated impulse responses to 1 percent productivity shock and 25 basis point monetary shock:
var ea = 1^2;    // unit shock to technology
var ev = 0.25^2; // shock of 25 basis points, i.e. 1 percentage point annualized
end;

//stoch_simul(irf=20);   
stoch_simul(irf=20, order=1) gap pi i r y n w v a;      



