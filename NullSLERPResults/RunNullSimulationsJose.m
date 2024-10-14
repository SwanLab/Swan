% RunNullSimulations

s.gJFlow = 0.7;
s.Vf     = 0.2;
Sim1     = ThreeDimCantileverDensity(s);

s.gJFlow = 2;
s.Vf     = 0.4;
Sim2     = ThreeDimCantileverDensity(s);