function Pol = polarization_ani(matprop)


ce = matprop.ce;
ci = matprop.ci;

[~,posdef]=chol(ce);

if posdef ~=0
    disp('constitutive tensor for matrix CE is not positive definited');
    stopprog;
end

[~,posdef]=chol(ci);

if posdef ~=0
    disp('constitutive tensor for inclusion CI is not positive definited');
    stopprog;
end

%% Polarization tensor for the bulk fase

invCe=inv(ce);
invCi=inv(ci);

pol = [invCe(1,1) -2*invCe(1,3) 2*invCe(1,2)+invCe(3,3) -2*invCe(2,3) invCe(2,2)];
Sol = roots(pol);

mu1 = Sol(1);
mu2 = Sol(3);

p1 = invCe(1,1).*mu1.^2 + invCe(1,2) - invCe(1,3)*mu1;
p2 = invCe(1,1).*mu2.^2 + invCe(1,2) - invCe(1,3)*mu2;
q1 = invCe(2,1).*mu1 + invCe(2,2)./mu1 - invCe(2,3);
q2 = invCe(2,1).*mu2 + invCe(2,2)./mu2 - invCe(2,3);

Q = (p1 - p2)./(2*(mu1 - mu2));
R = (p1*mu2 - mu1*p2)./(2*(mu1-mu2)); 
P = (q1 - q2)./(2*(mu1 - mu2));
O = (q1*mu2 - mu1*q2)./(2*(mu1-mu2)); 

RR = real(R);
RI = imag(R);
QR = real(Q);
QI = imag(Q);
PR = real(P);
PI = imag(P);
OR = real(O);
OI = imag(O);

CI11=invCi(1,1);
CI12=invCi(1,2);
CI13=invCi(1,3);
CI21=invCi(2,1);
CI22=invCi(2,2);
CI23=invCi(2,3);
CI31=invCi(3,1);
CI32=invCi(3,2);
CI33=invCi(3,3);

P11 = (2*(-(CI23^2*QI) + CI22*CI33*QI + 2*CI22*OI*QI + 2*CI33*OI*QI + 4*OI^2*QI + 4*OR^2*QI - 2*CI23*PI*QI - 4*OR*PI*QI + 2*CI22*PR*QI + 4*OI*PR*QI + 2*CI22*QI^2 + 4*OI*QI^2 - 2*CI22*PI*QR - 4*OI*PI*QR - 2*CI23*PR*QR - 4*OR*PR*QR + 2*CI22*QR^2 + 4*OI*QR^2 - 2*CI22*PI*RI - 4*OI*PI*RI - 2*CI23*PR*RI - 4*OR*PR*RI - 2*CI23*QI*RI + 4*OR*QI*RI - 4*PI*QI*RI + 2*CI22*QR*RI + 4*OI*QR*RI - 4*PR*QR*RI - 4*PR*RI^2 - CI13*((CI22 + 2*OI)*(PI - QR) + PR*(CI23 + 2*(OR + RI))) + CI12*((CI23 - 2*OR + 2*PI)*(PI - QR) + PR*(CI33 + 2*(OI + PR + QI - RR))) + 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 4*PR*RR^2))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P12 = (-2*(CI13*(-(CI23*OI) + CI22*(OR + RI)) + 2*(QR + RI)*(-(CI23*OI) + CI22*(OR + RI)) + CI12*(-((CI23 - 2*OR + 2*PI)*(OR + RI)) + OI*(CI33 + 2*(OI + PR + QI - RR))) + (CI23*(CI23 - 2*OR + 2*PI) - CI22*(CI33 + 2*(OI + PR + QI)))*RR + 2*CI22*RR^2))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P13 = (2*((QR + RI)*(-(CI33*(CI22 + 2*OI)) + CI23*(CI23 + 2*(OR + RI))) - CI12*(CI33*(OR - PI) + CI23*(OI + PR + QI - RR)) + CI13*((OR - PI)*(CI23 + 2*(OR + RI)) + (CI22 + 2*OI)*(OI + PR + QI - RR)) - 2*(CI33*(OR - PI) + CI23*(OI + PR + QI))*RR + 2*CI23*RR^2))/(-(CI13^2*(CI22 + 2*OI)) - CI12^2*(CI33 + 2*OI) + CI11*(-CI23^2 + 4*(OR - PI)*(OR + RI) - 2*CI23*(PI + RI) + (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*CI13*(-(CI22*(PI + RI)) + CI12*(CI23 + PI + RI) + CI23*(-PR + RR) - 2*OR*(PR + RR) - 2*(PR*RI + OI*(PI + RI) - PI*RR)) + 2*(-(CI23^2*QI) + CI22*CI33*QI + 2*CI22*OI*QI + 2*CI33*OI*QI + 4*OI^2*QI + 4*OR^2*QI - 4*OR*PI*QI + 2*CI22*PR*QI + 4*OI*PR*QI + 2*CI22*QI^2 + 4*OI*QI^2 - 2*CI22*PI*QR - 4*OI*PI*QR - 4*OR*PR*QR + 2*CI22*QR^2 + 4*OI*QR^2 - 2*CI22*PI*RI - 4*OI*PI*RI - 4*OR*PR*RI + 4*OR*QI*RI - 4*PI*QI*RI + 2*CI22*QR*RI + 4*OI*QR*RI - 4*PR*QR*RI - 4*PR*RI^2 + CI12*(CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI)) + (CI33 + 2*(OI + PR + QI - RR))*(PR - RR)) - CI12^2*(PR + QI - RR) + 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI - 2*(OR - PI)*(PI - QR))*RR - 4*PR*RR^2 - 2*CI23*(QI*RI + PR*(QR + RI) + PI*(QI - RR) + QR*RR)));
P21 = (-2*(-(CI13*((CI23 - 2*OR + 2*PI)*QI + PR*(CI13 + 2*(QR + RI)))) + CI11*((CI23 - 2*OR + 2*PI)*(PI - QR) + PR*(CI33 + 2*(OI + PR + QI - RR))) + CI12*(-((PI - QR)*(CI13 + 2*(QR + RI))) + QI*(CI33 + 2*(OI + PR + QI - RR)))))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P22 = (2*(-(CI13^2*OI) + 2*(CI33*OI*QI + 2*OI^2*QI + (OR + RI)*(-((CI23 - 2*OR + 2*PI)*QI) + (CI12 - 2*PR)*(QR + RI)) + 2*OI*(QI*(PR + QI) - (PI - QR)*(QR + RI))) + CI11*(CI33*OI - (CI23 - 2*OR + 2*PI)*(OR + RI) + 2*OI*(OI + PR + QI - RR)) - ((CI12 - 2*PR)*(CI33 + 2*(OI + PR)) + 2*(CI12 + 2*OI - 2*PR)*QI - 2*(CI23 - 2*OR + 2*PI)*(PI - QR))*RR + 2*(CI12 - 2*PR)*RR^2 + CI13*(CI12*(OR + RI) - 2*PR*(OR + RI) - 2*OI*(PI + RI) + (CI23 - 2*OR + 2*PI)*RR)))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P23 = (-2*(-2*CI13*OR*PI + 2*CI13*PI^2 + CI13^2*(-OR + PI) + 2*CI13*OI*PR + 2*CI13*PR^2 + 2*CI23*OI*QI + 2*CI33*OR*QI - 2*CI33*PI*QI + 2*CI13*PR*QI + 2*CI23*PR*QI + 2*CI23*QI^2 - CI13*CI23*QR + 2*CI13*OR*QR - 2*CI13*PI*QR - 2*CI23*PI*QR - 2*CI33*PR*QR + 2*CI23*QR^2 - CI13*CI23*RI - 2*CI23*PI*RI - 2*CI33*PR*RI + 2*CI23*QR*RI + CI12*(CI33*(QR + RI) - CI13*(OI + PR + QI - RR)) + CI11*(CI33*(OR - PI) + CI23*(OI + PR + QI - RR)) - 2*(CI13*PR + CI23*QI)*RR))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P31 = (2*(CI12^2*(-PI + QR) + CI11*((CI22 + 2*OI)*(PI - QR) + PR*(CI23 + 2*(OR + RI))) - CI13*(CI22*QI + 2*OI*QI + 2*PR*RR) + CI12*(-(CI13*PR) + QI*(CI23 + 2*(OR + RI)) + 2*(-PI + QR)*RR)))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P32 = (2*(-(CI12^2*(OR + RI)) - (CI11 + 2*QI)*(CI23*OI - CI22*(OR + RI)) - (2*CI23*PR + CI22*(CI13 + 2*PI - 2*QR))*RR + CI12*(OI*(CI13 + 2*PI - 2*QR) + 2*PR*(OR + RI) + CI23*RR)))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P33 = (2*(-2*CI22*OI*QI - 4*OI^2*QI - 2*CI23*OR*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + CI13*CI22*QR + 2*CI13*OI*QR + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + CI13*CI22*RI + 2*CI13*OI*RI + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI11*((CI22 + 2*OI)*(OI + PR + QI) + (OR - PI)*(CI23 + 2*(OR + RI))) + CI12*(CI13*(OR - PI) + 2*OR*(PI - 2*QR - RI) - CI23*(QR + RI) - 2*(PI^2 - PI*QR + RI*(QR + RI)) - 2*(PR - RR)*(OI + PR + QI - RR)) + CI12^2*(OI + PR + QI - RR) + CI11*(CI22 + 2*OI)*RR + 2*(-2*PR*(OI + PR) + (CI22 + 2*OI - 2*PR)*QI + (OR - PI)*(CI13 + 2*PI - 2*QR))*RR + 4*PR*RR^2))/(-(CI13^2*(CI22 + 2*OI)) - CI12^2*(CI33 + 2*OI) + CI11*(-CI23^2 + 4*(OR - PI)*(OR + RI) - 2*CI23*(PI + RI) + (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*CI13*(-(CI22*(PI + RI)) + CI12*(CI23 + PI + RI) + CI23*(-PR + RR) - 2*OR*(PR + RR) - 2*(PR*RI + OI*(PI + RI) - PI*RR)) + 2*(-(CI23^2*QI) + CI22*CI33*QI + 2*CI22*OI*QI + 2*CI33*OI*QI + 4*OI^2*QI + 4*OR^2*QI - 4*OR*PI*QI + 2*CI22*PR*QI + 4*OI*PR*QI + 2*CI22*QI^2 + 4*OI*QI^2 - 2*CI22*PI*QR - 4*OI*PI*QR - 4*OR*PR*QR + 2*CI22*QR^2 + 4*OI*QR^2 - 2*CI22*PI*RI - 4*OI*PI*RI - 4*OR*PR*RI + 4*OR*QI*RI - 4*PI*QI*RI + 2*CI22*QR*RI + 4*OI*QR*RI - 4*PR*QR*RI - 4*PR*RI^2 + CI12*(CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI)) + (CI33 + 2*(OI + PR + QI - RR))*(PR - RR)) - CI12^2*(PR + QI - RR) + 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI - 2*(OR - PI)*(PI - QR))*RR - 4*PR*RR^2 - 2*CI23*(QI*RI + PR*(QR + RI) + PI*(QI - RR) + QR*RR)));

PP = [P11 P12 P13; P21 P22 P23; P31 P32 P33];

deltaC = ci - ce;
Pe = deltaC*invCi*((eye(3) + PP)*ci*invCe - PP);

Pol.Pe = Pe;
Pol.PP = PP;
%% Polarization tensor for the inclusion

invCe=inv(ci);
invCi=inv(ce);

pol = [invCe(1,1) -2*invCe(1,3) 2*invCe(1,2)+invCe(3,3) -2*invCe(2,3) invCe(2,2)];
Sol = roots(pol);

mu1 = Sol(1);
mu2 = Sol(3);

p1 = invCe(1,1).*mu1.^2 + invCe(1,2) - invCe(1,3)*mu1;
p2 = invCe(1,1).*mu2.^2 + invCe(1,2) - invCe(1,3)*mu2;
q1 = invCe(2,1).*mu1 + invCe(2,2)./mu1 - invCe(2,3);
q2 = invCe(2,1).*mu2 + invCe(2,2)./mu2 - invCe(2,3);

Q = (p1 - p2)./(2*(mu1 - mu2));
R = (p1*mu2 - mu1*p2)./(2*(mu1-mu2)); 
P = (q1 - q2)./(2*(mu1 - mu2));
O = (q1*mu2 - mu1*q2)./(2*(mu1-mu2)); 

RR = real(R);
RI = imag(R);
QR = real(Q);
QI = imag(Q);
PR = real(P);
PI = imag(P);
OR = real(O);
OI = imag(O);

CI11=invCi(1,1);
CI12=invCi(1,2);
CI13=invCi(1,3);
CI21=invCi(2,1);
CI22=invCi(2,2);
CI23=invCi(2,3);
CI31=invCi(3,1);
CI32=invCi(3,2);
CI33=invCi(3,3);

P11 = (2*(-(CI23^2*QI) + CI22*CI33*QI + 2*CI22*OI*QI + 2*CI33*OI*QI + 4*OI^2*QI + 4*OR^2*QI - 2*CI23*PI*QI - 4*OR*PI*QI + 2*CI22*PR*QI + 4*OI*PR*QI + 2*CI22*QI^2 + 4*OI*QI^2 - 2*CI22*PI*QR - 4*OI*PI*QR - 2*CI23*PR*QR - 4*OR*PR*QR + 2*CI22*QR^2 + 4*OI*QR^2 - 2*CI22*PI*RI - 4*OI*PI*RI - 2*CI23*PR*RI - 4*OR*PR*RI - 2*CI23*QI*RI + 4*OR*QI*RI - 4*PI*QI*RI + 2*CI22*QR*RI + 4*OI*QR*RI - 4*PR*QR*RI - 4*PR*RI^2 - CI13*((CI22 + 2*OI)*(PI - QR) + PR*(CI23 + 2*(OR + RI))) + CI12*((CI23 - 2*OR + 2*PI)*(PI - QR) + PR*(CI33 + 2*(OI + PR + QI - RR))) + 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 4*PR*RR^2))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P12 = (-2*(CI13*(-(CI23*OI) + CI22*(OR + RI)) + 2*(QR + RI)*(-(CI23*OI) + CI22*(OR + RI)) + CI12*(-((CI23 - 2*OR + 2*PI)*(OR + RI)) + OI*(CI33 + 2*(OI + PR + QI - RR))) + (CI23*(CI23 - 2*OR + 2*PI) - CI22*(CI33 + 2*(OI + PR + QI)))*RR + 2*CI22*RR^2))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P13 = (2*((QR + RI)*(-(CI33*(CI22 + 2*OI)) + CI23*(CI23 + 2*(OR + RI))) - CI12*(CI33*(OR - PI) + CI23*(OI + PR + QI - RR)) + CI13*((OR - PI)*(CI23 + 2*(OR + RI)) + (CI22 + 2*OI)*(OI + PR + QI - RR)) - 2*(CI33*(OR - PI) + CI23*(OI + PR + QI))*RR + 2*CI23*RR^2))/(-(CI13^2*(CI22 + 2*OI)) - CI12^2*(CI33 + 2*OI) + CI11*(-CI23^2 + 4*(OR - PI)*(OR + RI) - 2*CI23*(PI + RI) + (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*CI13*(-(CI22*(PI + RI)) + CI12*(CI23 + PI + RI) + CI23*(-PR + RR) - 2*OR*(PR + RR) - 2*(PR*RI + OI*(PI + RI) - PI*RR)) + 2*(-(CI23^2*QI) + CI22*CI33*QI + 2*CI22*OI*QI + 2*CI33*OI*QI + 4*OI^2*QI + 4*OR^2*QI - 4*OR*PI*QI + 2*CI22*PR*QI + 4*OI*PR*QI + 2*CI22*QI^2 + 4*OI*QI^2 - 2*CI22*PI*QR - 4*OI*PI*QR - 4*OR*PR*QR + 2*CI22*QR^2 + 4*OI*QR^2 - 2*CI22*PI*RI - 4*OI*PI*RI - 4*OR*PR*RI + 4*OR*QI*RI - 4*PI*QI*RI + 2*CI22*QR*RI + 4*OI*QR*RI - 4*PR*QR*RI - 4*PR*RI^2 + CI12*(CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI)) + (CI33 + 2*(OI + PR + QI - RR))*(PR - RR)) - CI12^2*(PR + QI - RR) + 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI - 2*(OR - PI)*(PI - QR))*RR - 4*PR*RR^2 - 2*CI23*(QI*RI + PR*(QR + RI) + PI*(QI - RR) + QR*RR)));
P21 = (-2*(-(CI13*((CI23 - 2*OR + 2*PI)*QI + PR*(CI13 + 2*(QR + RI)))) + CI11*((CI23 - 2*OR + 2*PI)*(PI - QR) + PR*(CI33 + 2*(OI + PR + QI - RR))) + CI12*(-((PI - QR)*(CI13 + 2*(QR + RI))) + QI*(CI33 + 2*(OI + PR + QI - RR)))))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P22 = (2*(-(CI13^2*OI) + 2*(CI33*OI*QI + 2*OI^2*QI + (OR + RI)*(-((CI23 - 2*OR + 2*PI)*QI) + (CI12 - 2*PR)*(QR + RI)) + 2*OI*(QI*(PR + QI) - (PI - QR)*(QR + RI))) + CI11*(CI33*OI - (CI23 - 2*OR + 2*PI)*(OR + RI) + 2*OI*(OI + PR + QI - RR)) - ((CI12 - 2*PR)*(CI33 + 2*(OI + PR)) + 2*(CI12 + 2*OI - 2*PR)*QI - 2*(CI23 - 2*OR + 2*PI)*(PI - QR))*RR + 2*(CI12 - 2*PR)*RR^2 + CI13*(CI12*(OR + RI) - 2*PR*(OR + RI) - 2*OI*(PI + RI) + (CI23 - 2*OR + 2*PI)*RR)))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P23 = (-2*(-2*CI13*OR*PI + 2*CI13*PI^2 + CI13^2*(-OR + PI) + 2*CI13*OI*PR + 2*CI13*PR^2 + 2*CI23*OI*QI + 2*CI33*OR*QI - 2*CI33*PI*QI + 2*CI13*PR*QI + 2*CI23*PR*QI + 2*CI23*QI^2 - CI13*CI23*QR + 2*CI13*OR*QR - 2*CI13*PI*QR - 2*CI23*PI*QR - 2*CI33*PR*QR + 2*CI23*QR^2 - CI13*CI23*RI - 2*CI23*PI*RI - 2*CI33*PR*RI + 2*CI23*QR*RI + CI12*(CI33*(QR + RI) - CI13*(OI + PR + QI - RR)) + CI11*(CI33*(OR - PI) + CI23*(OI + PR + QI - RR)) - 2*(CI13*PR + CI23*QI)*RR))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P31 = (2*(CI12^2*(-PI + QR) + CI11*((CI22 + 2*OI)*(PI - QR) + PR*(CI23 + 2*(OR + RI))) - CI13*(CI22*QI + 2*OI*QI + 2*PR*RR) + CI12*(-(CI13*PR) + QI*(CI23 + 2*(OR + RI)) + 2*(-PI + QR)*RR)))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P32 = (2*(-(CI12^2*(OR + RI)) - (CI11 + 2*QI)*(CI23*OI - CI22*(OR + RI)) - (2*CI23*PR + CI22*(CI13 + 2*PI - 2*QR))*RR + CI12*(OI*(CI13 + 2*PI - 2*QR) + 2*PR*(OR + RI) + CI23*RR)))/(CI13^2*(CI22 + 2*OI) + CI12^2*(CI33 + 2*OI) + CI11*(CI23^2 - 4*(OR - PI)*(OR + RI) + 2*CI23*(PI + RI) - (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*(CI23^2*QI - CI22*CI33*QI - 2*CI22*OI*QI - 2*CI33*OI*QI - 4*OI^2*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI + 2*CI23*QI*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI12*(PR*(CI33 + 2*(OI + PR + QI)) + CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI))) + CI12^2*(PR + QI - RR) + CI12*(CI33 + 2*(OI + 2*PR + QI))*RR - 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI + (CI23 - 2*OR + 2*PI)*(PI - QR))*RR - 2*CI12*RR^2 + 4*PR*RR^2) + 2*CI13*(CI22*(PI + RI) - CI12*(CI23 + PI + RI) + CI23*(PR - RR) + 2*OR*(PR + RR) + 2*(PR*RI + OI*(PI + RI) - PI*RR)));
P33 = (2*(-2*CI22*OI*QI - 4*OI^2*QI - 2*CI23*OR*QI - 4*OR^2*QI + 2*CI23*PI*QI + 4*OR*PI*QI - 2*CI22*PR*QI - 4*OI*PR*QI - 2*CI22*QI^2 - 4*OI*QI^2 + CI13*CI22*QR + 2*CI13*OI*QR + 2*CI22*PI*QR + 4*OI*PI*QR + 2*CI23*PR*QR + 4*OR*PR*QR - 2*CI22*QR^2 - 4*OI*QR^2 + CI13*CI22*RI + 2*CI13*OI*RI + 2*CI22*PI*RI + 4*OI*PI*RI + 2*CI23*PR*RI + 4*OR*PR*RI - 4*OR*QI*RI + 4*PI*QI*RI - 2*CI22*QR*RI - 4*OI*QR*RI + 4*PR*QR*RI + 4*PR*RI^2 - CI11*((CI22 + 2*OI)*(OI + PR + QI) + (OR - PI)*(CI23 + 2*(OR + RI))) + CI12*(CI13*(OR - PI) + 2*OR*(PI - 2*QR - RI) - CI23*(QR + RI) - 2*(PI^2 - PI*QR + RI*(QR + RI)) - 2*(PR - RR)*(OI + PR + QI - RR)) + CI12^2*(OI + PR + QI - RR) + CI11*(CI22 + 2*OI)*RR + 2*(-2*PR*(OI + PR) + (CI22 + 2*OI - 2*PR)*QI + (OR - PI)*(CI13 + 2*PI - 2*QR))*RR + 4*PR*RR^2))/(-(CI13^2*(CI22 + 2*OI)) - CI12^2*(CI33 + 2*OI) + CI11*(-CI23^2 + 4*(OR - PI)*(OR + RI) - 2*CI23*(PI + RI) + (CI22 + 2*OI)*(CI33 + 2*(OI + PR + QI - RR))) + 2*CI13*(-(CI22*(PI + RI)) + CI12*(CI23 + PI + RI) + CI23*(-PR + RR) - 2*OR*(PR + RR) - 2*(PR*RI + OI*(PI + RI) - PI*RR)) + 2*(-(CI23^2*QI) + CI22*CI33*QI + 2*CI22*OI*QI + 2*CI33*OI*QI + 4*OI^2*QI + 4*OR^2*QI - 4*OR*PI*QI + 2*CI22*PR*QI + 4*OI*PR*QI + 2*CI22*QI^2 + 4*OI*QI^2 - 2*CI22*PI*QR - 4*OI*PI*QR - 4*OR*PR*QR + 2*CI22*QR^2 + 4*OI*QR^2 - 2*CI22*PI*RI - 4*OI*PI*RI - 4*OR*PR*RI + 4*OR*QI*RI - 4*PI*QI*RI + 2*CI22*QR*RI + 4*OI*QR*RI - 4*PR*QR*RI - 4*PR*RI^2 + CI12*(CI23*(PI + RI) + 2*(PI^2 - PI*QR + RI*(QR + RI) + OR*(-PI + 2*QR + RI)) + (CI33 + 2*(OI + PR + QI - RR))*(PR - RR)) - CI12^2*(PR + QI - RR) + 2*(PR*(CI33 + 2*(OI + PR)) - (CI22 + 2*OI - 2*PR)*QI - 2*(OR - PI)*(PI - QR))*RR - 4*PR*RR^2 - 2*CI23*(QI*RI + PR*(QR + RI) + PI*(QI - RR) + QR*RR)));

PP = [P11 P12 P13; P21 P22 P23; P31 P32 P33];

deltaC = ce - ci;
Pi = deltaC*invCi*((eye(3) + PP)*ce*invCe - PP);

Pol.Pi = Pi;
end

