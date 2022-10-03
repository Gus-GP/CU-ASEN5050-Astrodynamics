%Gustavo Grinsteins
%ASEN 5050
%HW6 - Defining Lamberts Equation

function F = LambertEqt(mu,a,s,c,TOF,less180,greatThanTOFmin)

alpha = 2*asin(sqrt((s)/(2*a)));
beta = 2*asin(sqrt((s-c)/(2*a)));
n = sqrt(mu/a^3);

%perform sign checks
if less180
    %beta = beta;
else
    beta = -beta;
end

if greatThanTOFmin
    alpha = 2*pi - alpha;
else
    %alpha = alpha;
end

%define the entire equation set to zero
F = (1/n)*((alpha-beta)-(sin(alpha)-sin(beta)))-TOF;

