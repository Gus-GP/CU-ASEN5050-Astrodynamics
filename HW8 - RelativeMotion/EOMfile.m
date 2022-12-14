function [statedot] = EOMfile(t,state,mu)
x = state(1);
y = state(2);
z = state(3);
R = [x,y,z];
r = norm(R);
xdot = state(4);
ydot = state(5);
zdot = state(6);
statedot(1) = xdot;
statedot(2) = ydot;
statedot(3) = zdot;
statedot(4) = -(mu/(r^3))*x;
statedot(5) = -(mu/(r^3))*y;
statedot(6) = -(mu/(r^3))*z;
statedot = statedot';