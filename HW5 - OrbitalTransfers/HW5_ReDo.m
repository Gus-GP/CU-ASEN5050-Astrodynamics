%House keeping
clc;
clear;
%Problem 1
% mu_moon = 4902.799; %km^3/s^2
% a = 8500;
% e = 0.29;
% theta_star = -21*(pi/180);
% delta_V = [0.25, -0.06, 0];
% 
% p = a*(1-e^2);
% h = sqrt(mu_moon*p);
% 
% v_r = (mu_moon/h)*e*sin(theta_star);
% v_theta = (mu_moon/h)*(1 + e*cos(theta_star));
% 
% V_1 = [v_r,v_theta,0];
% 
% V_2 = V_1 + delta_V;
% v_2 = norm(V_2);
% 
% r_1 = p/(1+e*cos(theta_star));
% r_2 = r_1;
% 
% h_2 = r_2*V_2(2);
% mech_e_2 = (v_2^2)/(2) - (mu_moon)/(r_2);
% a_2 = -(mu_moon)/(mech_e_2*2);
% e_2 = sqrt(1 + (2*mech_e_2*h_2^2)/(mu_moon^2));
% p_2 = a_2*(1-e_2^2);
% theta_star_2 = abs(acos((p_2-r_2)/(r_2*e_2)));
% if dot([r_2,0,0],V_2) < 0
%     theta_star_2 = -theta_star_2;
% end
% 
% delta_omega = theta_star - theta_star_2;

%% Problem 2
% mu_mars = 4.305*10^4; %km^3/s^2
% EQR_mars = 3397.2; %km
% mech_e_bef_1 = -5.16157; %km^2/s^2
% E_bef_1 = -1.46057; %rad
% ts_bef_1 = -90*(pi/180); %rad
% 
% a_bef_1 = -(mu_mars)/(2*mech_e_bef_1);
% ecc_bef_1 = cos(E_bef_1);
% p_bef_1 = a_bef_1*(1-ecc_bef_1^2);
% h_bef_1 = sqrt(mu_mars*p_bef_1);
% 
% v_r_bef_1 = (mu_mars/h_bef_1)*ecc_bef_1*sin(ts_bef_1);
% v_theta_bef_1 = (mu_mars/h_bef_1)*(1+ecc_bef_1*cos(ts_bef_1));
% 
% V_bef_1 = [v_r_bef_1,v_theta_bef_1,0];
% 
% rp_bef_1 = a_bef_1*(1-ecc_bef_1);
% altp_bef_1 = rp_bef_1 - EQR_mars;
% 
% ra_bef_1 = a_bef_1*(1+ecc_bef_1);
% 
% r_1 = p_bef_1/(1+ecc_bef_1*cos(ts_bef_1));
% 
% %same apo diff peri
% rp_after_1 = 400 + EQR_mars;
% ra_after_1 = ra_bef_1;
% a_aft_1 = (rp_after_1 + ra_after_1)/2;
% mech_e_aft_1 = - (mu_mars)/(2*a_aft_1);
% ecc_aft_1 = 1 - (rp_after_1/a_aft_1);
% p_aft_1 = a_aft_1*(1-ecc_aft_1^2);
% h_aft_1 = sqrt(mu_mars*p_aft_1);
% ts_aft_1 = -abs(acos((p_aft_1-r_1)/(r_1*ecc_aft_1)));
% 
% v_r_aft_1 = (mu_mars/h_aft_1)*ecc_aft_1*sin(ts_aft_1);
% v_theta_aft_1 = (mu_mars/h_aft_1)*(1+ecc_aft_1*cos(ts_aft_1));
% 
% V_aft_1 = [v_r_aft_1,v_theta_aft_1,0];
% 
% maneuver = V_aft_1 - V_bef_1;
% 
% %given maneuver
% deltaV = [0.04,-0.002,0];
% 
% V2 = V_bef_1 + deltaV;
% 
% mech_e_new = (norm(V2)^2)/2 - mu_mars/r_1;
% 
% a_new = -(mu_mars)/(2*mech_e_new);
% 
% h_new = r_1*V2(2);
% 
% ecc_new = sqrt(1 + (2*mech_e_new*h_new^2)/(mu_mars^2));
% 
% %new altitude
% AltitudeA_new = (a_new*(1-ecc_new)) - EQR_mars

%% Problem 3
mu_sun = 1.32712428*10^11; %km^3/s^2
a_earth = 1.0000010178*149597870.7;
a_jupi = 9.554909595*149597870.7;

at = (a_earth+a_jupi)/2;

v_i = sqrt(mu_sun/a_earth);
v_pt = sqrt(((2*mu_sun)/(a_earth)) - (mu_sun/at));

deltaV1 = v_pt - v_i;%in the theta hat direction

v_f = sqrt(mu_sun/a_jupi);
v_at = sqrt(((2*mu_sun)/(a_jupi)) - (mu_sun/at));

deltaV2 = v_f - v_at;

total_deltav = deltaV1 + deltaV2;
TOF = pi*sqrt((at^3)/mu_sun);

%calculate phase angle to randezvous
n_jupi = sqrt(mu_sun/(a_jupi^3));
alpha = n_jupi*TOF;
phi = pi - alpha;

%part c
rB = 11*149597870.7;

ab_t_1 = (rB + a_earth)/2;

vb_pt_1 = sqrt(((2*mu_sun)/(a_earth)) - (mu_sun/ab_t_1));

deltaVb_1 = vb_pt_1 - v_i;

vb_at_1 = sqrt(((2*mu_sun)/(rB)) - (mu_sun/ab_t_1));

ab_t_2 = (rB + a_jupi)/2;

vb_at_2 = sqrt(((2*mu_sun)/(rB)) - (mu_sun/ab_t_2));

deltaVb_2 = vb_at_2 - vb_at_1;

vb_pt_2 = sqrt(((2*mu_sun)/(a_jupi)) - (mu_sun/ab_t_2));

deltaVb_3 = v_f - vb_pt_2;

total_deltavb = abs(deltaVb_1) + abs(deltaVb_2) + abs(deltaVb_3);

TOFb = pi*sqrt((ab_t_1^3)/mu_sun) + pi*sqrt((ab_t_2^3)/mu_sun);
