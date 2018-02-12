function heavisidePerimeter
syms epsilon;
syms x;
syms txi
syms long;
long_value = 1;
a = -long_value;
b = long_value;


alpha_e = 0.5 + 1/pi*atan(txi);
alpha_e2 = 1/(1+exp(-2*txi));
%alpha_e_amstutz = -0.5*(exp(2*txi_max)-1)/(exp(4*txi_max)-1)*exp(txi) + 0.5*(exp(2*txi_max)-1)/(exp(4*txi_max)-1)*exp(-txi);

epsil = [1,1e-1,1e-2,1e-3];
for ieps = 1:length(epsil)
   alph_e_atan = matlabFunction(alpha_e);
   alpha_e_exp = matlabFunction(alpha_e2);
    
   txi_max_plot = long_value/epsil(ieps);
   alpha_e_ams = @(txi) alpha_e_amstutz(txi,txi_max_plot);
   txi_plot = linspace(-txi_max_plot,txi_max_plot,100);
   xplot = txi_plot*epsil(ieps);
   plot(xplot,alph_e_atan(txi_plot),xplot,alpha_e_exp(txi_plot),xplot,alpha_e_ams(txi_plot))
    
end


Per_epsilon = 0.5*subs(int(1-alpha_e),txi,long_value/epsilon);
Per_epsilon = matlabFunction(Per_epsilon);

Per_epsilon2 = 0.5*subs(int(1-alpha_e2),txi,long_value/epsilon);
Per_epsilon2 = matlabFunction(Per_epsilon2);


eps = 10.^[1:-1:-5];

for ieps = 1:length(eps)
    epsil = eps(ieps);
    Per_epsilon_value(ieps) = Per_epsilon(epsil);
    Per_epsilon_value2(ieps) = Per_epsilon2(epsil);
    
    txi_max_plot = long_value/epsil;
    alpha_e_ams = @(txi) 1-alpha_e_amstutz(txi,txi_max_plot);
    Per_e_amstutz(ieps) = 0.5*(integral(alpha_e_ams,0,long_value/epsil));
end

semilogx(eps,Per_epsilon_value,eps,Per_epsilon_value2,eps,Per_e_amstutz)


end


function alpha_e_amstutz = alpha_e_amstutz(txi,txi_max)
alpha_e_amstutz1 = zeros(size(txi));   
alpha_e_amstutz1(txi < 0) = 0.5*exp(txi(txi < 0));
alpha_e_amstutz1(txi > 0) = 1-0.5*exp(-txi(txi > 0));

alpha_e_amstutz2 = -0.5*(exp(2*txi_max)-1)/(exp(4*txi_max)-1)*exp(txi) + 0.5*(exp(2*txi_max)-1)/(exp(4*txi_max)-1)*exp(-txi);
alpha_e_amstutz = alpha_e_amstutz1 ;%+ alpha_e_amstutz2;
end