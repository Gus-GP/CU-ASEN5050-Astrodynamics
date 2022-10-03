%Gustavo Grinsteins
%ASEN 5050
%HW2 Problem 2

%House Keeping
clc;
clear;

%% Part A
fprintf('Problem 2 Part A \n')
%Position Vector
R = [4981.75,-4121.90,22.70]; %Km
%Velocity Vector 
V = [-0.60359,0.56812,-2.24093]; %Km/s
%Gravitational Constant
mu = 4.305*10^4; %Km^3/s^2

%magnitudes
r = norm(R); %Km
v = norm(V); %Km/s
fprintf('r is %4.2f Km, v is %1.4f Km/s \n\n',r,v)

%Calculating specific angular momentum
H = cross(R,V); %Km^2/s
h = norm(H); %Km^2/s
fprintf('H in XYZ frame is <%4.2f,%2.2f,%4.2f> Km^2/s\n',H)
fprintf('h is %4.2f Km^2/s \n\n',h)

%Calculating Specific Energy
Sp_E = ((v)^2)/(2) - (mu/r); %Km^2/s^2
fprintf('Specific Energy = %4.4f Km^2/s^2 \n\n',Sp_E)


%Inclination angle
Zhat = [0,0,1];
i = max(min(dot(H,Zhat)/(norm(H)*norm(Zhat)),1),-1);%Radians
iDegrees = real(acosd(i));%Degrees
fprintf('inclination angle i = %4.4f deg \n\n',iDegrees)

%Semi-major axis
a = -mu/(2*Sp_E);%Km
fprintf('semi-major axis a = %4.2f Km\n\n',a)

%Eccentricity Vector
Ecc = cross(V,H)*(1/mu) - R/r;%Unitless
ecc = norm(Ecc);%Unitless

fprintf('Eccentricity vector in XYZ frame is <%4.4f,%2.4f,%4.4f> \n',Ecc)
fprintf('e is %4.4f\n\n',ecc)

%RAAN
N = cross(Zhat,H);%Km^2/s
n = norm(N);%Km^2/s
fprintf('N vector in XYZ frame is <%4.4f,%4.4f,%4.4f> Km^2/s\n',N)
fprintf('n is %4.4f Km^2/s\n\n',n)
Xhat = [1,0,0];
Om = max(min(dot(N,Xhat)/(norm(N)*norm(Xhat)),1),-1);%Radians
OmDegrees = real(acosd(Om));%Degrees
fprintf('RAAN angle Omega = %4.4f deg \n\n',OmDegrees)

%Argument of periapsis
w = max(min(dot(N,Ecc)/(norm(N)*norm(Ecc)),1),-1);%Radians
wDegrees = -1*real(acosd(w));%Degrees
fprintf('Argument of Periapsis w = %4.4f deg \n\n',wDegrees)

%True Anomaly
ThetaStar = max(min(dot(R,Ecc)/(norm(R)*norm(Ecc)),1),-1);%Radians
TSDegrees = -1*real(acosd(ThetaStar));%Degrees
fprintf('True Anomaly ThetaStar = %4.4f deg \n\n',TSDegrees)

%Rotation matrix
theta = TSDegrees+wDegrees;%Degrees

%Rotation matrix
R1 = [1,0,0;0,cosd(iDegrees),sind(iDegrees);0,-sind(iDegrees),cosd(iDegrees)];
R3_Om = [cosd(OmDegrees),sind(OmDegrees),0;-sind(OmDegrees),cosd(OmDegrees),0;0,0,1];
R3_theta = [cosd(theta),sind(theta),0;-sind(theta),cosd(theta),0;0,0,1];

C = R3_theta*R1*R3_Om;

%% Part B
fprintf('Problem 2 Part B \n')
%Transforming position from XYZ to r,theta,h
Rrot = C*transpose(R);
fprintf('R vector in (r,theta,h) frame is <%4.4f,%4.4f,%4.4f> Km\n\n',Rrot)

%Transforming velocity from XYZ to r,theta,h
Vrot = C*transpose(V);
fprintf('V vector in (r,theta,h) frame is <%4.4f,%4.4f,%4.4f> Km/s\n\n',Vrot)

%% Part C S/C at ascending node
fprintf('Problem 2 Part C \n')
PositionRot = [3904.4447;0;0];

VelocityRot = [0.8368;3.7073;0];

%Update rotation matrix
%Rotation matrix
theta = 0;%Degrees

%Rotation matrix
R1 = [1,0,0;0,cosd(iDegrees),sind(iDegrees);0,-sind(iDegrees),cosd(iDegrees)];
R3_Om = [cosd(OmDegrees),sind(OmDegrees),0;-sind(OmDegrees),cosd(OmDegrees),0;0,0,1];
R3_theta = [cosd(theta),sind(theta),0;-sind(theta),cosd(theta),0;0,0,1];

C = R3_theta*R1*R3_Om;

%Transforming position from r,theta,h to XYZ 
PositionXYZ = C.'*PositionRot;
fprintf('R vector in XYZ frame is <%4.4f,%4.4f,%4.4f> Km\n\n',PositionXYZ)

%Transforming velocity from r,theta,h to XYZ
VelocityXYZ = C.'*VelocityRot;
fprintf('V vector in XYZ frame is <%4.4f,%4.4f,%4.4f> Km/s\n\n',VelocityXYZ)




