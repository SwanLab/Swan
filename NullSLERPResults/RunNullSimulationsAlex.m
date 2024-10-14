% RunNullSimulations

s.gJFlow = 0.5;
s.Vf     = 0.4;
Sim1     = ThreeDimCantileverDensity(s);

s.gJFlow = 1;
s.Vf     = 0.4;
Sim2     = ThreeDimCantileverDensity(s);