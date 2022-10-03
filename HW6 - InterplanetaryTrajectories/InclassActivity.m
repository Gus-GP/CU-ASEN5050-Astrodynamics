%Gustavo Grinsteins
%ASEN 5050
%HW6

%House Keeping
clc;
clear;
%% In class activity check
%Given
mu_earth = 3.986004415*10^5;
%earth vectors
R1 = [-654,13605,1997];
V1 = [-5.53,0.849,0.6830];
%earth vectors
R2 = [7284,-19341,-3264];
V2 = [3.07,2.63,0.444];
r1 = norm(R1);
r2 = norm(R2);
%STEP 1: Find the transfer angle (theta>180)
less180 = false;
delta_ThetaStar = abs(acos(dot(R1,R2)/(r1*r2)));
delta_ThetaStar = 2*pi-delta_ThetaStar;
%STEP 2: Calculate Geometric quantities
c = sqrt(r1^2+r2^2-2*(r1*r2)*cos(delta_ThetaStar));
s = (1/2)*(r1+r2+c);
%STEP 3: Determine if the TOF is for an ellipse
TOFp = (1/3)*sqrt(2/mu_earth)*((s^(3/2))+((s-c)^(3/2)));
TOF = 5*60*60;
if TOF > TOFp
    fprintf('The transfer orbit is elliptical \n')
else
    fprintf('The transfer orbit is not elliptical \n')
end
%STEP 4: Determine correct alpha and beta
a_m = (s/2);
n_m = sqrt(mu_earth/(a_m^3));
alpha_m = pi;
beta_m_0 = 2*asin(sqrt((s-c)/(s)));
if less180
    beta_m = beta_m_0;%(theta<180)
else
    beta_m = -beta_m_0;%(theta>180)
end
TOFmin = (1/n_m)*((alpha_m-beta_m)-(sin(alpha_m)-sin(beta_m)));

a = EllOrbitLambertEqSolve(s,c,a_m,TOF,TOFmin,mu_earth,less180);

n = sqrt(mu_earth/(a^3));
alpha_0 = 2*asin(sqrt((s)/(2*a)));
beta_0 = 2*asin(sqrt((s-c)/(2*a)));
if less180
    beta = beta_0;
else
    beta = -beta_0;
end
if TOF > TOFmin
    alpha = 2*pi - alpha_0;
else
    alpha = alpha_0;
end
TOF_new = (1/n)*((alpha-beta)-(sin(alpha)-sin(beta)));