function rho_eng = rho_eng_fun(M,rho,Earth_cons)
%Constants:
gamma = Earth_cons(1,8); %Gamma of air

%Calculations:
Base = 1+0.5.*(gamma-1).*M.^2;
exponent = -1./(gamma-1);

Denominator = Base.^exponent;

rho_eng = rho./Denominator;
end

