function minimize_coef_of_Trace_top_opt

problem.objective =  @(x) beta_function(x);
problem.lb = [-1,-1];
problem.ub = [1/2,1/2];
problem.options = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point','MaxIter',1000000,'MaxFunEvals',1000000,...
                                    'TolFun',1e-12,'TolX',1e-12);
problem.solver = 'fmincon';
problem.x0 = [0.3,0.3];

fmincon(problem)

end


function beta = beta_function(x)
Ee = 1;
Ei = 10;
nue = x(1);
nui = x(2);

beta = (- 3*Ee^2*nue*nui^2 + 3*Ee^2*nue + Ee^2*nui^2 - Ee^2 + 6*Ee*Ei*nue^2*nui - 2*Ee*Ei*nue*nui + 2*Ee*Ei*nue - 8*Ee*Ei*nui + 2*Ee*Ei - 3*Ei^2*nue^3 + Ei^2*nue^2 + 3*Ei^2*nue - Ei^2)/((nue^2 - 1)*(- Ee^2*nui^2 + Ee^2 + 2*Ee*Ei*nue*nui - 2*Ee*Ei*nui + 4*Ee*Ei - Ei^2*nue^2 + 2*Ei^2*nue + 3*Ei^2));


end