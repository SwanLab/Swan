function CL_alpha = CL_alpha_fun(M,Aircraft_cons,Aero_cons)
%Constants:
Cl_alpha_0 = Aero_cons(1,2);
e_factor = Aero_cons(1,3);
AR = Aircraft_cons(1,3);

%Calculations:

Denominator = 1 + (Cl_alpha_0)./(e_factor.*AR.*pi);

CL_alpha_0 = Cl_alpha_0./Denominator;

%Calculation with variable CL_alpha(M)

Denominator = sqrt(1-M^2);

CL_alpha = CL_alpha_0./Denominator;
end

