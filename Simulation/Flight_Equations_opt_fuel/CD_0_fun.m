function CD_0 = CD_0_fun(M,Aero_cons)
%Constants
a = Aero_cons(1,4);
b = Aero_cons(1,5);
c = Aero_cons(1,6);
d = Aero_cons(1,7);
e = Aero_cons(1,8);

%Calculations:
CD_0 = a + b.*M + c.*M.^2 + d.*M.^3 + e.*M.^4;
end

