%Gustavo Grinsteins
%ASEN 5050
%HW4 - Iterative function

function ComputedE = NewtonRaphsonMethodForE(e,n,elapsedTimeSec)
%NewtonRaphsonMethod: To solve Kepler's equation

%initial values
M = n*elapsedTimeSec;%Radians
E0 = M; %Initial value for eccentric anomaly
diff = 1000;%Initial value for tolerance stopping condition
tolerance = 10^-4; %4 digit accuracy is desired
epsilon = 10^-14; %Do not divide by a number smaller than this

%Iterate until the tolerance codition is met
    while (abs(diff) > tolerance)
    f = M - E0 + e*sin(E0);
    fprime = - 1 + e*cos(E0);
    %divide by zero check
        if (abs(fprime)< epsilon)
            fprintf('divide by zero exception \n\n')
            %If f prime is getting to small break from the loop
            break
        end
    %New eccentric anomaly calc by newton-raphson method
    E1 = E0 - (f)/(fprime);
    %update values to keep iterating
    diff = abs(E1-E0);
    E0 = E1;
    end 
ComputedE = E1;
fprintf('The eccentric anomaly after %4.2f minutes pass periapsis is %4.4f radians\n',elapsedTimeSec*(1/60),ComputedE)
end

