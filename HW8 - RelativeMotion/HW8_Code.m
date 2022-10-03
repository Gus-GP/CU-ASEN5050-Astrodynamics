%Gustavo Grinsteins
%ASEN 5050
%HW8

%House Keeping
clc;
clear;

%% Problem 1
%given
mu_sun = 1.32712428*10^11;
G = 6.673*10^-20;
rp = 7500;
ra = 8500;
i = 105*(pi/180);
Period = 110*(60);%seconds
rP = 6500;
aP = 2.25*149597870.7;
a = (rp+ra)/2;
ecc = 1 - (rp/a);
mu_P = (a^3)*(((2*pi)/(Period))^2);
mP = mu_P/G;
vP = sqrt(mu_sun/aP);
PlanetOrbitTime = 2*pi*aP*(1/vP);
RAAN_dot = (2*pi)/(PlanetOrbitTime);
J2 = -(2/3)*(RAAN_dot/cos(i))*(((1-ecc^2)^2*a^(7/2))/(sqrt(mu_P)*rP^2));

%%  Problem 2
%given
mu_earth = 3.986004415*10^5;
a_earth = 1.0000010178*149597870.7;
R0 = [-6402,-1809,1065];
V0 = [0.999,-6.471,-4.302];

%magnitudes
r0 = norm(R0); %Km
v0 = norm(V0); %Km/s
H = cross(R0,V0); %Km^2/s
h = norm(H); %Km^2/s
p = (h^2)/mu_earth;
mech_e = ((v0)^2)/(2) - (mu_earth/r0); %Km^2/s^2
%Semi-major axis
a = -mu_earth/(2*mech_e);%Km
%Eccentricity Vector
Ecc = cross(V0,H)*(1/mu_earth) - R0/r0;%Unitless
ecc = norm(Ecc);%Unitless
theta_star0 = abs(acos((p-r0)/(r0*ecc)));
if dot(V0,R0) < 0
   theta_star0 = -theta_star0; 
end
%calculatime time from tp to t0
E_0 = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(theta_star0/2));
n = sqrt(mu_earth/(a^3));
t0_minus_tp = (1/n)*(E_0 - ecc*sin(E_0));
tp_to_t1 = t0_minus_tp + (60*60);
E_1 = NewtonRaphsonMethodForE(ecc,n,tp_to_t1);
r1 = a*(1-ecc*cos(E_1));
%calulate state vector with f and g func
deltaE = E_1-E_0;
f_1 = 1 - (a/r0)*(1-cos(deltaE));%unitless
g_1 = (tp_to_t1-t0_minus_tp) - sqrt((a^3)/mu_earth)*(deltaE-sin(deltaE));%seconds
fdot_1 = (-sin(deltaE)*sqrt(mu_earth*a))/(r0*r1);%1/s
gdot_1 = 1 - (a/r1)*(1-cos(deltaE));%s
R1 = f_1*R0 + g_1*V0;
V1 = fdot_1*R0 + gdot_1*V0;
h1 = norm(cross(R1,V1));
mech_e_1 = ((norm(V1))^2)/(2) - (mu_earth/norm(R1)); %Km^2/s^2
tp_to_t2 = t0_minus_tp + (100*60*60);
E_2 = NewtonRaphsonMethodForE(ecc,n,tp_to_t2);
r2 = a*(1-ecc*cos(E_2));
deltaE = E_2-E_0;
f_2 = 1 - (a/r0)*(1-cos(deltaE));%unitless
g_2 = (tp_to_t2-t0_minus_tp) - sqrt((a^3)/mu_earth)*(deltaE-sin(deltaE));%seconds
fdot_2 = (-sin(deltaE)*sqrt(mu_earth*a))/(r0*r2);%1/s
gdot_2 = 1 - (a/r2)*(1-cos(deltaE));%s
R2 = f_2*R0 + g_2*V0;
V2 = fdot_2*R0 + gdot_2*V0;
h2 = norm(cross(R2,V2));
mech_e_2 = ((norm(V2))^2)/(2) - (mu_earth/norm(R2)); %Km^2/s^2
%% ODE stuff
state0 = [-6402,-1809,1065,0.999,-6.471,-4.302];
%1
options = odeset('Stats','off','RelTol',1*10^-12,'AbsTol',1*10^-12);
[tout1,xout1] = ode45(@EOMfile,[0 60*60],state0,options,mu_earth);
output1 = xout1(end,:);
R1ode = [output1(1),output1(2),output1(3)];
V1ode = [output1(4),output1(5),output1(6)];
partCdeltaR1 = norm(R1ode-R1);
partCdeltaV1 = norm(V1ode-V1);
h1ode = norm(cross(R1ode,V1ode));
mech_e_1ode = ((norm(V1ode))^2)/(2) - (mu_earth/norm(R1ode)); %Km^2/s^2
%2
options = odeset('Stats','off','RelTol',1*10^-12,'AbsTol',1*10^-12);
[tout2,xout2] = ode45(@EOMfile,[0 100*60*60],state0,options,mu_earth);
output2 = xout2(end,:);
output2beg = xout2(1,:);
R2ode12beg = [output2beg(1),output2beg(2),output2beg(3)];
V2ode12beg = [output2beg(4),output2beg(5),output2beg(6)];
h2ode12beg = norm(cross(R2ode12beg,V2ode12beg));
mech_e_2ode12beg = ((norm(V2ode12beg))^2)/(2) - (mu_earth/norm(R2ode12beg));
R2ode12 = [output2(1),output2(2),output2(3)];
V2ode12 = [output2(4),output2(5),output2(6)];
partCdeltaR2 = norm(R2ode12-R2);
partCdeltaV2 = norm(V2ode12-V2);
h2ode12 = norm(cross(R2ode12,V2ode12));
mech_e_2ode12 = ((norm(V2ode12))^2)/(2) - (mu_earth/norm(R2ode12)); %Km^2/s^2
deltaR_12 = norm(R2ode12-R2);
deltaV_12 = norm(V2ode12-V2);
deltah_12 = h2ode12-h2ode12beg;
deltaMech_e_12 = mech_e_2ode12-mech_e_2ode12beg;
%%%%%%
options = odeset('Stats','off','RelTol',1*10^-10,'AbsTol',1*10^-10);
[tout2,xout2] = ode45(@EOMfile,[0 100*60*60],state0,options,mu_earth);
output2 = xout2(end,:);
output2beg = xout2(1,:);
R2ode10beg = [output2beg(1),output2beg(2),output2beg(3)];
V2ode10beg = [output2beg(4),output2beg(5),output2beg(6)];
h2ode10beg = norm(cross(R2ode10beg,V2ode10beg));
mech_e_2ode10beg = ((norm(V2ode10beg))^2)/(2) - (mu_earth/norm(R2ode10beg));
R2ode10 = [output2(1),output2(2),output2(3)];
V2ode10 = [output2(4),output2(5),output2(6)];
h2ode10 = norm(cross(R2ode10,V2ode10));
mech_e_2ode10 = ((norm(V2ode10))^2)/(2) - (mu_earth/norm(R2ode10)); %Km^2/s^2
deltaR_10 = norm(R2ode10-R2);
deltaV_10 = norm(V2ode10-V2);
deltah_10 = h2ode10-h2ode10beg;
deltaMech_e_10 = mech_e_2ode10-mech_e_2ode10beg;
%%%%%%%
options = odeset('Stats','off','RelTol',1*10^-8,'AbsTol',1*10^-8);
[tout2,xout2] = ode45(@EOMfile,[0 100*60*60],state0,options,mu_earth);
output2 = xout2(end,:);
output2beg = xout2(1,:);
R2ode8beg = [output2beg(1),output2beg(2),output2beg(3)];
V2ode8beg = [output2beg(4),output2beg(5),output2beg(6)];
h2ode8beg = norm(cross(R2ode8beg,V2ode8beg));
mech_e_2ode8beg = ((norm(V2ode8beg))^2)/(2) - (mu_earth/norm(R2ode8beg));
R2ode8 = [output2(1),output2(2),output2(3)];
V2ode8 = [output2(4),output2(5),output2(6)];
h2ode8 = norm(cross(R2ode8,V2ode8));
mech_e_2ode8 = ((norm(V2ode8))^2)/(2) - (mu_earth/norm(R2ode8)); %Km^2/s^2
deltaR_8 = norm(R2ode8-R2);
deltaV_8 = norm(V2ode8-V2);
deltah_8 = h2ode8-h2ode8beg;
deltaMech_e_8 = mech_e_2ode8-mech_e_2ode8beg;
%%%%%%%
options = odeset('Stats','off','RelTol',1*10^-6,'AbsTol',1*10^-6);
[tout2,xout2] = ode45(@EOMfile,[0 100*60*60],state0,options,mu_earth);
output2 = xout2(end,:);
output2beg = xout2(1,:);
R2ode6beg = [output2beg(1),output2beg(2),output2beg(3)];
V2ode6beg = [output2beg(4),output2beg(5),output2beg(6)];
h2ode6beg = norm(cross(R2ode6beg,V2ode6beg));
mech_e_2ode6beg = ((norm(V2ode6beg))^2)/(2) - (mu_earth/norm(R2ode6beg));
R2ode6 = [output2(1),output2(2),output2(3)];
V2ode6 = [output2(4),output2(5),output2(6)];
h2ode6 = norm(cross(R2ode6,V2ode6));
mech_e_2ode6 = ((norm(V2ode6))^2)/(2) - (mu_earth/norm(R2ode6)); %Km^2/s^2
deltaR_6 = norm(R2ode6-R2);
deltaV_6 = norm(V2ode6-V2);
deltah_6 = h2ode6-h2ode6beg;
deltaMech_e_6 = mech_e_2ode6-mech_e_2ode6beg;
%%%%%%%
options = odeset('Stats','off','RelTol',1*10^-4,'AbsTol',1*10^-4);
[tout2,xout2] = ode45(@EOMfile,[0 100*60*60],state0,options,mu_earth);
output2 = xout2(end,:);
output2beg = xout2(1,:);
R2ode4beg = [output2beg(1),output2beg(2),output2beg(3)];
V2ode4beg = [output2beg(4),output2beg(5),output2beg(6)];
h2ode4beg = norm(cross(R2ode4beg,V2ode4beg));
mech_e_2ode4beg = ((norm(V2ode4beg))^2)/(2) - (mu_earth/norm(R2ode4beg));
R2ode4 = [output2(1),output2(2),output2(3)];
V2ode4 = [output2(4),output2(5),output2(6)];
h2ode4 = norm(cross(R2ode4,V2ode4));
mech_e_2ode4 = ((norm(V2ode4))^2)/(2) - (mu_earth/norm(R2ode4)); %Km^2/s^2
deltaR_4 = norm(R2ode4-R2);
deltaV_4 = norm(V2ode4-V2);
deltah_4 = h2ode4-h2ode4beg;
deltaMech_e_4 = mech_e_2ode4-mech_e_2ode4beg;