%Gustavo Grinsteins
%ASEN 5050
%HW9

%House Keeping
clc;
clear;

%problem 1
mu_earth = 3.986004415*10^5;
r = 10000;
n = sqrt(mu_earth/r^3);
for t = 0:10000
    x1(t+1) = -15*cos(n*t);
    y1(t+1) = 15 + 30*sin(n*t);
end
figure(1)
plot(y1,x1)
%Problem 2
TOF = (2/3)*pi*sqrt(r^3/mu_earth);
t_0 = 0;
t_1 = TOF;
R0 = [-5;5;0];
V_0_minus = [0.04;-0.01;0.01];
R1 = [5;5;0];
V_1_plus = [0;0;0];
%defining CW matrix eqs
phi_rr_t1 = [4-3*cos(n*t_1),0,0;6*(sin(n*t_1)-n*t_1),1,0;0,0,cos(n*t_1)];
phi_rv_t1 = [(1/n)*sin(n*t_1),(2/n)*(1-cos(n*t_1)),0;(2/n)*(cos(n*t_1)-1),(4/n)*sin(n*t_1)-3*t_1,0;0,0,(1/n)*sin(n*t_1)];
phi_vr_t1 = [3*n*sin(n*t_1),0,0;6*(n*cos(n*t_1)-n),0,0;0,0,-n*sin(n*t_1)];
phi_vv_t1 = [cos(n*t_1),2*sin(n*t_1),0;-2*sin(n*t_1),4*cos(n*t_1)-3,0;0,0,cos(n*t_1)];

V_0_plus = -inv(phi_rv_t1)*(phi_rr_t1*R0-R1);
deltaV1 = V_0_plus-V_0_minus;
deltav1 = norm(deltaV1);

V_1_minus = phi_vr_t1*R0+phi_vv_t1*V_0_plus;
deltaV2 = V_1_plus-V_1_minus;
deltav2 = norm(deltaV2);

R0_2 = [-5;5];
V_0_minus_2 = [0.04;-0.01];
V_0_plus_2 = [V_0_plus(1);V_0_plus(2)];

x_0 = -5;
y_0 = 5;
xdot_0 = 0.04;
ydot_0 = -0.01;

for t_i = 0:TOF
    placeholder = [4-3*cos(n*t_i),0;6*(sin(n*t_i)-n*t_i),1]*R0_2+[(1/n)*sin(n*t_i),(2/n)*(1-cos(n*t_i));(2/n)*(cos(n*t_i)-1),(4/n)*sin(n*t_i)-3*t_i]*V_0_plus_2;
    %x(t_i+1) = 4*x_0+(2/n)*(ydot_0)+(xdot_0/n)*sin(n*t_i)-(3*x_0+(2/n)*ydot_0)*cos(n*t_i);
    %y(t_i+1) = y_0-(2/n)*xdot_0-3*(2*n*x_0+ydot_0)*t_i+2*(3*x_0+(2/n)*ydot_0)*sin(n*t_i)+(2/n)*xdot_0*cos(n*t_i);
    x(t_i+1) = placeholder(1);
    y(t_i+1) = placeholder(2);
end
figure(2)
plot(y,x)

%% Problem 2
mu_sun = 1.32712428*10^11;
mu_mars = 4.305*10^4;
mars_EQR = 3397.2;
J2_mars = 0.001964;
a_sc = mars_EQR + 300;
a_mars = 1.52367934*149597870.7;
v_mars = sqrt(mu_sun/a_mars);
circOrbTime = (2*pi*a_mars)*(1/v_mars);
%period = 2*pi*sqrt(a_mars^3/mu_sun);
RAAN_dot = (2*pi)/circOrbTime;

incli = acos((-RAAN_dot*(2/3)*(a_sc^(7/2)))/(sqrt(mu_mars)*J2_mars*mars_EQR^2));

%b
aero_period = 1.02595675*86400;
a_sc_aero = ((aero_period/(2*pi))^(2)*mu_mars)^(1/3);
altitude_b = a_sc_aero - mars_EQR;
RAAN_dot_aero = -((3/2)*((sqrt(mu_mars)*J2_mars*mars_EQR^2)/(a_sc_aero^(7/2))));



