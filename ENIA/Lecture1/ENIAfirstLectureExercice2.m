function ENIAfirstLectureExercice2
g = 9.81;
l = 1;
w = sqrt(g/l);
m = 1;

theta0V  = pi/20;
dtheta0V = 0;

tV = linspace(0,5,1000);


thetaV_E = obtainAnaliticalSolutionViaExplicitFormula(tV,w,theta0V,dtheta0V);
plot(tV,thetaV_E,'+')

thetaV_A = obtainAnaliticalSolutionViaSymbolicProcedure(tV,w,theta0V,dtheta0V);
hold on
plot(tV,thetaV_A,'+')

error = norm(thetaV_E - thetaV_A)/norm(thetaV_A)

dthetaV = obtainDerivativeViaExplicitFormula(tV,w,theta0V,dtheta0V);
figure
plot(tV,dthetaV,'+')

[E,Ek,Ep] = computeEnergies(m,g,l,thetaV_A,dthetaV);
figure(3)
plot(tV,[E;Ek;Ep])

end

function thetaV = obtainAnaliticalSolutionViaExplicitFormula(tV,w,theta0V,dtheta0V)
B = dtheta0V/w;
A = theta0V;
thetaV = A*cos(w*tV) + B*sin(w*tV);
end

function dthetaV = obtainDerivativeViaExplicitFormula(tV,w,theta0V,dtheta0V)
B = dtheta0V/w;
A = theta0V;
dthetaV = w*A*sin(w*tV) - w*B*cos(w*tV);
end

function [E,Ek,Ep] = computeEnergies(m,g,l,theta,dtheta)
Ek = computeKineticEnergy(m,dtheta);
Ep = computePotentialEnergy(m,g,l,theta);
E = Ek + Ep;
end

function Ep = computePotentialEnergy(m,g,l,theta)
Ep = m*g*l*(1-cos(theta));
end

function Ek = computeKineticEnergy(m,dtheta)
Ek = 0.5*m*dtheta.*dtheta;
end


function thetaV = obtainAnaliticalSolutionViaSymbolicProcedure(tV,w,theta0V,dtheta0V)
A = sym('A');
B = sym('B');
%
t = sym('t');
%
theta0 = sym('theta0','real');
dtheta0 = sym('dtheta0','real');
%
theta = A*cos(w*t) + B*sin(w*t);
%
dtheta = diff(theta,t);
%
eq(1) = theta0  == subs(theta,t,0);
eq(2) = dtheta0 == subs(dtheta,t,0);
%
solution = solve(eq,[A,B]);
%
theta = subs(theta,[A,B],[solution.A,solution.B]);
theta = simplify(theta);
%

theta = subs(theta,theta0,theta0V);
theta = subs(theta,dtheta0,dtheta0V);

theta = matlabFunction(theta);

thetaV = theta(tV);

end