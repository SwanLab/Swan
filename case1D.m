% 1D case for Continuum Damage Model 
%% Linear hardening law
clc,clear,close all

E0 = 210;
H  = -0.5;
r0 = 10;
r1 = 20;
qInf = r0 + H*(r1-r0);
rOld = r0;

uVals = [0:1e-1:10];
converged = false;

qVec = [];
F = [];
dmg = []; 
tauVec = [];

for i=1:length(uVals)
    u = uVals(i);
    %while(~converged)
        % Compute R
        tau = sqrt(u*E0*u);
        if tau>rOld
            r = tau;
        else
            r = rOld;
        end

        % Compute d
        if r>=r1
            q = qInf;
        else
            q = r0 + H*(r-r0);
        end

        d = 1 - q/r;
        if d>1
            d = 1;
        elseif d<0
            d = 0;
        end

        % Compute reactions
        reac = (1-d)*E0*u;
    %end
    rOld = r;

    dmg(end+1) = d;
    qVec(end+1) = q;
    F(end+1) = reac;
    tauVec(end+1) = tau;
end

figure()
plot(tauVec,dmg)

figure()
plot(tauVec,qVec)

figure()
plot(tauVec,F)

%% Exponential hardening law
clc,clear,close all

E0 = 210;
A = 0.1;
r0 = 10;
r1 = 20;
H  = 0.5;
qInf = -(r0 + H*(r1-r0));
rOld = r0;

uVals = [0:1e-1:10];
converged = false;

qVec = [];
F = [];
dmg = []; 
tauVec = [];

for i=1:length(uVals)
    u = uVals(i);
    %while(~converged)
        % Compute R
        tau = sqrt(u*E0*u);
        if tau>rOld
            r = tau;
        else
            r = rOld;
        end

        % Compute d
        q = qInf - (qInf - r0)*exp(A*(1-r/r0));

        d = 1 - q/r;
        if d>1
            d = 1;
        elseif d<0
            d = 0;
        end

        % Compute reactions
        reac = (1-d)*E0*u;
    %end
    rOld = r;

    dmg(end+1) = d;
    qVec(end+1) = q;
    F(end+1) = reac;
    tauVec(end+1) = tau;
end

figure()
plot(tauVec,dmg)

figure()
plot(tauVec,qVec)

figure()
plot(tauVec,F)