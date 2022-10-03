%Gustavo Grinsteins
%ASEN 5050
%HW5

%House Keeping
clc;
clear;
%% Problem 1
%Given
mu_moon = 4902.799;
EQR_moon = 1738;
a1 = 8500;
e1 = 0.29;
theta_star_1 = -21*(pi/180);
p1 = a1*(1-e1^2);
h1 = sqrt(mu_moon*p1);
vr1 = (mu_moon/h1)*e1*sin(theta_star_1);
vtheta1 = (mu_moon/h1)*(1+e1*cos(theta_star_1));
vrthetah1 = [vr1,vtheta1,0];
vdelta = [0.25,-0.06,0];
vrthetah2 = vrthetah1 + vdelta;
r1 = (p1)/(1+e1*cos(theta_star_1));
mech_e2 = (((norm(vrthetah2))^2)/2) - (mu_moon)/(r1);
h2 = r1*vrthetah2(2);
e2 = sqrt(1+((2*mech_e2*h2^2)/(mu_moon^2)));
a2 = -(mu_moon)/(2*mech_e2);
p2 = a2*(1-e2^2);
theta_star_2c = -acos((p2-r1)/(r1*e2));
theta_star_2 = -acos((vrthetah2(2)*h2 - mu_moon)/(mu_moon*e2));
%% Problem 2
%Given
mu_mars = 4.305*10^4;
mars_EQR = 3397.2;
mech_e_bef = -5.16187;
E1_bef = -1.46057;
thetas_1bef = -90*(pi/180);
e_bef = cos(E1_bef);
a_bef = -(mu_mars)/(2*mech_e_bef);
p_bef = a_bef*(1-e_bef^2);
h_bef = sqrt(mu_mars*p_bef);
vr1 = (mu_mars/h_bef)*e_bef*sin(thetas_1bef);
vtheta1 = (mu_mars/h_bef)*(1+e_bef*cos(thetas_1bef));
ra = a_bef*(1+e_bef);
rp_aft = mars_EQR + 400;
a_aft = (ra+rp_aft)/2;
e_aft = (ra-a_aft)/(a_aft);
p_aft = a_aft*(1-e_aft^2);
h_aft = sqrt(mu_mars*p_aft);
r_aft = p_bef;
theta_star_aft = -acos((p_aft-p_bef)/(p_bef*e_aft));
vr2 = (mu_mars/h_aft)*e_aft*sin(theta_star_aft);
vtheta2 = (mu_mars/h_aft)*(1+e_aft*cos(theta_star_aft));
deltaV = [vr2,vtheta2,0] - [vr1,vtheta1,0];
Vdelta_c = [0.04,-0.002,0];
V2rot = [vr1,vtheta1,0] + Vdelta_c;
mech_e_aft = (((norm(V2rot))^2)/2) - (mu_mars)/(p_bef);
a_aft_c = -(mu_mars/(2*mech_e_aft));
h_aft_c = p_bef*V2rot(2);
e_aft_c = sqrt(1-((h_aft_c ^2)/(mu_mars*a_aft_c)));
rp_aft_c = a_aft_c*(1-e_aft_c);
altitude_c = rp_aft_c - mars_EQR;
%% Problem 3
% %Given
mu_sun = 1.32712428*10^11;
AU_convert = 149597870.7;
a_i = 1.0000010178*AU_convert;%this is a for earth
a_f = 9.554909595*AU_convert;%This is a for saturn
r_pt = a_i; %transfer orbit periapsis radius
r_at = a_f; %transfer orbit apoapsis radius
v1_i = sqrt(mu_sun/a_i); %velocity at earths orbit
v2_f = sqrt(mu_sun/a_f); %velocity at saturns orbit
a_t = (1/2)*(a_i+a_f);%semi major axis of transfer orbit
e_t = (a_f-a_i)/(a_f+a_i);%eccentricity of transfer orbit
mech_e_t = -(mu_sun)/(2*a_t);%mech e of transfer orbit
vp_t = sqrt(((2*mu_sun)/(a_i))-(mu_sun/a_t));
va_t = sqrt(((2*mu_sun)/(a_f))-(mu_sun/a_t));
deltaV_1 = vp_t - v1_i;
deltaV_2 = v2_f - va_t;
totalDeltaV = deltaV_1+deltaV_2;
TOF = pi*sqrt((a_t^3)/mu_sun);
n_saturn = sqrt(mu_sun/(a_f^3));
alpha_l = TOF*n_saturn;
phase_angle = pi-alpha_l;
%PartC
rB = 11*AU_convert;
at1 = (1/2)*(a_i+rB);
et1 = (rB-a_i)/(rB+a_i);
at2 = (1/2)*(a_f+rB);
et2 = (rB-a_f)/(rB+a_f);
v1_i = sqrt(mu_sun/a_i);
vp_t1 = sqrt(((2*mu_sun)/(a_i))-(mu_sun/at1));
delta_V_1 = vp_t1 - v1_i;
va_t1 = sqrt(((2*mu_sun)/(rB))-(mu_sun/at1));
va_t2 = sqrt(((2*mu_sun)/(rB))-(mu_sun/at2));
delta_V_2 = va_t2 - va_t1;
vp_t2 = sqrt(((2*mu_sun)/(a_f))-(mu_sun/at2));
v3_f = sqrt(mu_sun/a_f);
delta_V_3 = v3_f - vp_t2;
totalDeltaV_2 = abs(delta_V_1)+abs(delta_V_2)+abs(delta_V_3);
TOF_bi = pi*sqrt((at1^3)/mu_sun)+pi*sqrt((at2^3)/mu_sun);
percentDiff = ((TOF_bi-TOF)/(TOF))*100;
WRH = 1 - exp((-totalDeltaV)/(301*.00981));
WRB = 1 - exp((-totalDeltaV_2)/(301*.00981));