function stress = StrainMultiplStress(strain,coeff,COMP) 
stress = coeff.*strain ; 
stress(COMP.c3) = 0.5*stress(COMP.c3) ; 


