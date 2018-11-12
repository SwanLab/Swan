function [P_mas,P_menos] = polarization_iso(gamma,mu,lambda)

[h1_menos,h2_menos] = coef_gmacro(gamma,mu,lambda);
[h1_mas,h2_mas] = coef_gmacro(1/gamma,mu,lambda);


P_mas = [h1_mas+h2_mas, h2_mas 0; h2_mas, h1_mas+h2_mas 0; 0 0 h1_mas];
P_menos = [h1_menos+h2_menos, h2_menos 0; h2_menos, h1_menos+h2_menos 0; 0 0 h1_menos];

end