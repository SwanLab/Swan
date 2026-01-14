clc
clear all 
syms z h tauXZ  eXZ  G eXZavg

eXZ = (1 - 4*z^2/h^2)*(3/2*eXZavg) ; 
tauXZ = G*eXZ ; 
Qx = int((1 - 4*z^2/h^2)*tauXZ,z,-h/2,h/2)

QxAVG = int(tauXZ,z,-h/2,h/2)