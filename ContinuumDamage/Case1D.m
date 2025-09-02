u=0:1e-3:1;
E=210;
r0=1/sqrt(6);
r1 = 5;
H=0,5;

Edamage = zeros(1,length(u));
Edamage(1) = r0;
sigma = zeros(1,length(u));
sigmaBar = zeros(1,length(u));

r = zeros(1,length(u));
r(1) = r0;
q = zeros(1,length(u));
d = zeros(1,length(u));
d_Dot = zeros(1,length(u));

qInf=r0 + H*(r1-r0);

for i =  2:length(u)

    epsi = u(i);
    sigBar= E*epsi;
    sigmaBar(i) = sigBar;
    tauEpsi = sqrt(epsi*sigBar);
    
    if ( tauEpsi > r(i-1))
        r(i) = tauEpsi;
    else
        r(i) = r0;
    end
   
    if ( r(i) <= r0 )
        d_Dot (i) = 0;
        q(i)= r0;
    elseif (r(i)>r0 && r(i) <r1)
        q(i) = r0 + H*(r(i)-r0);
        d_Dot(i) =(q(i)-H*r(i))/r(i)^3;
    else
        q(i) = qInf;
        d_Dot(i) =(q(i))/r(i)^3;
    end

d(i) = 1-q(i)/r(i);

Edamage(i) = (1-d(i))*E - d_Dot(i)*sigBar^2;
sig(i) = Edamage(i)*epsi;

end
close all
plot(r,d)

% Vectoritzar aixo quan es tingui temps!!!