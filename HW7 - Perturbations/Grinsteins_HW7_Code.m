%Gustavo Grinsteins
%ASEN 5050
%HW6

%House Keeping
clc;
clear;

%% Problem 1
%given
mu_sun = 1.32712428*10^11; %km^3/s^2
mu_earth = 3.986004415*10^5;
mu_mercury = 2.2032*10^4;
AU_to_km = 149597870.7; %km/AU
a_Earth = 1.0000010178*AU_to_km; %km
a_mercury = 0.387098309*AU_to_km; %km
r_mercury = 2439; %km
r_periapsis = 4000; %km

%Calcs
a_t = (1/2)*(a_Earth+a_mercury);
ecc = (a_Earth/a_t)-1;
Period = pi*sqrt(((a_t)^3)/(mu_sun));
R_SOI_Merc = ((mu_mercury/mu_sun)^(2/5))*a_mercury;
v_Merc = sqrt(mu_sun/a_mercury);
v_in = sqrt(((2*mu_sun)/a_mercury)-(mu_sun/a_t));
v_inf_in = v_in - v_Merc;
mech_e_hyp = v_inf_in^2/2;
a_hyp = -(mu_mercury)/(2*mech_e_hyp);
e_hyp = 1-(r_periapsis/a_hyp);
turning_angle = 2*asin(1/e_hyp);
v_out = sqrt(v_Merc^2+v_inf_in^2-2*v_Merc*v_inf_in*cos(pi-turning_angle));

%% Problem 2
mu_saturn = 3.794*10^7;
ra = 2000000;
rp = 700000;
a_titan = 1221830;
a_saturn = (1/2)*(ra+rp);
ecc = (ra/a_saturn)-1;
p = a_saturn*(1-ecc^2);
h = sqrt(mu_saturn*p);
theta_star_intersect = acos((p-a_titan)/(a_titan*ecc));
theta_star_b = acos(-ecc);
%Rotating frame at this theta star
vtheta = (mu_saturn/h)*(1+ecc*cos(theta_star_intersect));
vr = (mu_saturn/h)*ecc*sin(theta_star_intersect);
v_titan = sqrt(mu_saturn/a_titan);
V_in = [vr,vtheta,0];
V_titan = [0,v_titan,0];
v_inf_in = V_in - V_titan;
%Part b
r_pb = 2800;
mu_titan = (6.673*10^-20)*(1.3455*10^23); 
mech_e_hypb = (norm(v_inf_in)^2)/2;
a_hypb = -(mu_titan)/(2*mech_e_hypb);
e_hypb = 1-(r_pb/a_hypb);
turning_angleb = 2*asin(1/e_hypb);
%part c
%Angle between theta hat and vinf
beta = acos((dot(v_inf_in,[0,1,0]))/(norm(v_inf_in)));
v_inf_out = [norm(v_inf_in)*sin(beta-turning_angleb),norm(v_inf_in)*cos(beta-turning_angleb),0];
V_out = V_titan + v_inf_out;
%Part d
%assuming instantaneous change
mech_e_aft = ((norm(V_out)^2)/2) - (mu_saturn/a_titan);
a_aft = (-mu_saturn)/(2*mech_e_aft);
h_aft = a_titan*V_out(2);
ecc_aft = sqrt(1+((2*h_aft^2*mech_e_aft)/(mu_saturn^2)));
p_after = a_aft*(1-ecc_aft^2); 
theta_star_after = abs(acos((p_after-a_titan)/(a_titan*ecc_aft)));
if dot(V_out,[1,0,0]) < 0 
    theta_star_after = -theta_star_after;
end
%Part e
delta_V_eq = V_out - V_in;

%% Problem 3
mu_jupiter = 1.268*10^8; %km^3/s^2
PeriapsisAlt = 9.8999*10^6 - 71492; 
mech_e_hyp = (10.27299)^2/2;
a_hyp = -(mu_jupiter)/(2*mech_e_hyp);
e_hyp = 1-((9.8999*10^6)/a_hyp);
diff = 4.29053*10^6 - 60268;