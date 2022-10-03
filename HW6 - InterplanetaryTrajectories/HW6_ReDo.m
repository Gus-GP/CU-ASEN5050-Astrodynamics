%House Keeping
clc;
clear;

%% Problem 1
% mu_sun = 1.32712428*10^11;
% earth_julianDate = 2459528.5;
% R1 = [1.00078*10^8,1.09250*10^8,-5.29404*10^3];
% V1 = [-22.46086,20.00474,-1.79921*10^-4];
% venus_julianDate = 2459640.5;
% R2 = [-1.05048*10^8,-2.37576*10^7,5.73539*10^6];
% V2 = [7.49784,-34.31464,-0.90369];
% 
% %figuring out transfer details
% delta_TS = abs(acos(dot(R1,R2)/(norm(R1)*norm(R2))));
% TOF = (venus_julianDate-earth_julianDate)*86400;
% %Calculating geometric properties
% %chord length
% c = sqrt(norm(R1)^2+norm(R2)^2-2*norm(R1)*norm(R2)*cos(delta_TS));
% %semi-perimeter
% s = (0.5)*(norm(R1)+norm(R2)+c);
% %figuring out alpha and beta
% a_min = s/2;
% alpha_min = pi;
% beta_min = 2*asin(sqrt((s-c)/s));
% beta_min = beta_min;%since delta_ts <180
% n_min = sqrt(mu_sun/((a_min)^3));
% TOF_min = (1/n_min)*((alpha_min-beta_min)-(sin(alpha_min)-sin(beta_min)));
% %check that this is actually an ellipse
% TOF_pless = (1/3)*sqrt(2/mu_sun)*(s^(3/2)-((s-c)^(3/2)));
% if TOF > TOF_pless
%     fprintf("is an ellipse \n");
%     a_transfer = EllOrbitLambertEqSolve(s,c,a_min,TOF,TOF_min,mu_sun,true);
% end
% %recalculating alpha and beta values
% alpha_0 = 2*asin(sqrt(s/(2*a_transfer)));
% beta_0 = 2*asin(sqrt((s-c)/(2*a_transfer)));
% beta = beta_0; %because delta_TS_star < 180
% if TOF > TOF_min
%     alpha = 2*pi-alpha_0;
% else
%     alpha = alpha_0;
% end
% %confirm TOF
% n = sqrt(mu_sun/(a_transfer^3));
% TOF_calc = (1/n)*((alpha-beta)-(sin(alpha)-sin(beta)));
% %finding eccentricity
% p = (((4*a_transfer*(s-norm(R1))*(s-norm(R2))))/(c^2))*sin((alpha+beta)/2)^2;
% e = sqrt(1-(p/a_transfer));
% mech_e = -(mu_sun)/(2*a_transfer);
% %Part d calculating maneuvers
% V1_i = V1;
% V2_f = V2;
% 
% %calculating TSs
% thetaS_1 = -abs(acos((1/e)*((p/norm(R1))-1)));
% thetaS_2 = -abs(acos((1/e)*((p/norm(R2))-1)));
% 
% %calculating transfer velocities
% f = 1- (norm(R2)/p)*(1-cos(delta_TS));
% g = (norm(R2)*norm(R1)*sin(delta_TS))/(sqrt(mu_sun*p));
% f_dot = sqrt(mu_sun/p)*tan(delta_TS/2)*(((1-cos(delta_TS))/(p))-(1/norm(R2))-(1/norm(R1)));
% g_dot = 1-(norm(R1)/p)*(1-cos(delta_TS));
% 
% V1_f = (R2-f*R1)/(g);
% 
% V2_i = f_dot*R1+g_dot*V1_f;
% 
% deltaV_1 = V1_f - V1_i;
% deltav1 = norm(deltaV_1);
% deltaV_2 = V2_f - V2_i;
% deltav2 = norm(deltaV_2);

%% class problem
mu_earth = 3.986004415*10^5;
R1 = [-654,13605,1997];
V1 = [-5.53,0.849,0.6830];
R2 = [7284,-19341,-3264];
V2 = [3.07,2.63,0.444];
TOF = 5*60*60;% > 180
delta_TS = abs(acos(dot(R1,R2)/(norm(R1)*norm(R2))));
delta_TS = 2*pi-delta_TS;
c = sqrt(norm(R1)^2+norm(R2)^2-2*norm(R1)*norm(R2)*cos(delta_TS));
s = (0.5)*(norm(R1)+norm(R2)+c);
a_min = s/2;
alpha_min = pi;
beta_min = 2*asin(sqrt((s-c)/s));
beta_min = -beta_min;%since delta_ts >180
n_min = sqrt(mu_earth/((a_min)^3));
TOF_min = (1/n_min)*((alpha_min-beta_min)-(sin(alpha_min)-sin(beta_min)));
TOF_pless = (1/3)*sqrt(2/mu_earth)*(s^(3/2)+((s-c)^(3/2)));
if TOF > TOF_pless
    fprintf("is an ellipse \n");
    a_transfer = EllOrbitLambertEqSolve(s,c,a_min,TOF,TOF_min,mu_earth,false);
end
%recalculating alpha and beta values
alpha_0 = 2*asin(sqrt(s/(2*a_transfer)));
beta_0 = 2*asin(sqrt((s-c)/(2*a_transfer)));
beta = -beta_0; %because delta_TS_star > 180
if TOF > TOF_min
    alpha = 2*pi-alpha_0;
else
    alpha = alpha_0;
end
%confirm TOF
n = sqrt(mu_earth/(a_transfer^3));
TOF_calc = (1/n)*((alpha-beta)-(sin(alpha)-sin(beta)));
%finding eccentricity
p = (((4*a_transfer*(s-norm(R1))*(s-norm(R2))))/(c^2))*sin((alpha+beta)/2)^2;
e = sqrt(1-(p/a_transfer));
mech_e = -(mu_earth)/(2*a_transfer);
%Part d calculating maneuvers
V1_i = V1;
V2_f = V2;

%calculating TSs
thetaS_1 = -abs(acos((1/e)*((p/norm(R1))-1)));
thetaS_2 = -abs(acos((1/e)*((p/norm(R2))-1)));

%calculating transfer velocities
f = 1- (norm(R2)/p)*(1-cos(delta_TS));
g = (norm(R2)*norm(R1)*sin(delta_TS))/(sqrt(mu_earth*p));
f_dot = sqrt(mu_earth/p)*tan(delta_TS/2)*(((1-cos(delta_TS))/(p))-(1/norm(R2))-(1/norm(R1)));
g_dot = 1-(norm(R1)/p)*(1-cos(delta_TS));

V1_f = (R2-f*R1)/(g);

V2_i = f_dot*R1+g_dot*V1_f;

deltaV_1 = V1_f - V1_i;
deltav1 = norm(deltaV_1);
deltaV_2 = V2_f - V2_i;
deltav2 = norm(deltaV_2);