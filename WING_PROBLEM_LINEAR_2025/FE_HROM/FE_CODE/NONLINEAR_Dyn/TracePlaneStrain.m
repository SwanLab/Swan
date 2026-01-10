function trEident = TracePlaneStrain(eNP1,COMP)

 % COMP.c1 = 1:nstrain:ngausT ; 
% COMP.c2 = 2:nstrain:ngausT ; 
% COMP.c4 = 4:nstrain:ngausT ; 

trE = eNP1(COMP.c1) + eNP1(COMP.c2) + eNP1(COMP.c4);
trEident = zeros(size(eNP1)) ; 
trEident(COMP.c1) =   trE ; 
trEident(COMP.c2) =   trE ; 
trEident(COMP.c4) =   trE ; 