function strain = StressMultiplStrain(stress,coeff,COMP) 
strain = coeff.*stress ; 
strain(COMP.c3) = 2*strain(COMP.c3) ; 


