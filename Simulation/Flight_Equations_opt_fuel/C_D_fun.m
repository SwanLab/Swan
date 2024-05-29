function C_D = C_D_fun(CL_alpha,CD_0,ALPHA_Ajusted,Aero_cons)
%Constants:
K = Aero_cons(1,1);

C_D = CD_0 + K.*(CL_alpha.^2).*ALPHA_Ajusted.^2;
end

