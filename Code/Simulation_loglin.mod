// Begin 'copy_of_simulate_multiple_sector.mod'

@#ifndef a_val
    @#define a_val= -5
@#endif

@#ifndef theta_val
    @#define theta_val= 6
@#endif

@#ifndef rho_val
    @#define rho_val= 0
@#endif


var y y1 y2 c g1 g2 tau pi1 pi pi2 psi1 psi2 i s_g1 Gbar U; 
varexo eps_g;
parameters n s1p sp alpha1 alpha2 beta phi rho1 rho2 rho_g s1g a gamma1 gamma2 A1 A2 theta Lambda1 Lambda2 chi1 chi2 rho delta;


n      = 0.65;
s1p    = 0.716;
sp     = 0.813;
alpha1 = 0.78;
alpha2 = 0.89;
beta   = 0.997;
phi    = 4;
rho1   = 0.65;
rho2   = 0.73;
s1g    = 0.24;
rho_g  = 0.68;
gamma1 = n^(-phi);
gamma2 = (1-n)^(-phi);
mu     = 1
rho   = @{rho_val};
theta = @{theta_val};
a     = @{a_val};
A1     = 1-(s1g*(1-sp)/sp)*mu/(theta-1);
A2     = 1-((1-s1g)*(1-sp)/sp)*mu/(theta-1);
Lambda1 = (1 - alpha1)*(1 - beta*alpha1)/(alpha1*A1);
Lambda2 = (1 - alpha2)*(1 - beta*alpha2)/(alpha2*A2);
chi1    = s1g * (1-sp) / n;
chi2    = (1-s1g) * (1-sp)/ (1-n);
delta = -ln((1-s1g)/(s1g)); 

// Your model equations
model;
    // Aggregate output
    y = n*y1 + (1 - n)*y2;
    
    // Output in sector 1
    n*y1 = - s1p*sp*(1 - s1p)*tau + s1p*sp*c + s1g*(1 - sp)*g1;
    
    // Output in sector 2
    (1 - n)*y2 = (1 - s1p)*sp*s1p*tau + (1 - s1p)*sp*c + (1 - s1g)*(1 - sp)*g2;
    
    // New Keynesian Phillips Curve sector 1
    pi1 = (alpha1+ (1-alpha1)/A1)*beta*pi1(+1) + (1 - alpha1)*(1 - beta*alpha1)/(alpha1*A1)*(psi1 + (1-A1)*g1);
    
    // New Keynesian Phillips Curve sector 2
    pi2 = (alpha2+ (1-alpha2)/A2)*beta*pi2(+1) + (1 - alpha2)*(1 - beta*alpha2)/(alpha2*A2)*(psi2 + (1-A2)*g2);
    
    // Real product wage sector 1
    psi1 = c + phi*y1 - (1 - s1p)*tau;
    
    // Evolution of terms of trade
    tau = tau(-1) + pi1 - pi2;
    
    // Real product wage sector 2
    psi2 = c + phi*y2 + s1p*tau;
    
    // Euler equation
    c = c(+1) - i + pi(+1);
    
    // Taylor rule
    i = 1.5 * pi;
    
    // Inflation condition
    pi = s1p*pi1 + (1 - s1p)*pi2;
    
    // AR(1) process for total government spending
    Gbar = rho_g * Gbar(-1) + eps_g;
    
    // Logistic function for government spending share
    s_g1 = (delta + a)*(1-s1g)*Gbar;
    g1 = (1+(1-s1g)*(delta + a))*Gbar;
    (1-s1g)^2/(1-s1g + s1g^(1+rho)/(1-s1g)^rho)*g2 = Gbar - (s1g)^2/(s1g+((1-s1g)^(1+rho)/s1g^rho)) * g1;
    
    // Log-linearized utility function using steady-state values
    U = -(s1p*((1+phi)*(y1^2) + theta * (1-chi1)* pi1^2 / Lambda1) + (1-s1p)*((1+phi)*y2^2 + theta * (1-chi2) * pi2^2 / Lambda2));
end;

initval;
    y    = 1;
    y1   = 0.5;
    y2   = 0.5;
    c    = 1;
    g1   = 0;
    g2   = 0;
    tau  = 1;
    pi1  = 0;
    pi2  = 0;
    psi1 = 0;
    psi2 = 0;
    i    = 0;
    pi   = 0;
    Gbar = 0;
    U    = 0;
end;

shocks;
    var eps_g;    stderr 0.133;
end;

stoch_simul(order=2, periods=100000, irf=20, noprint);

// End 'copy_of_simulate_multiple_sector.mod'
