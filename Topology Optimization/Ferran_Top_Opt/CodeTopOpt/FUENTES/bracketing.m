phi = 30;
theta = 90;
Vol = 0.6;

lambda0 = 0.01;
vol_0 = MAIN_OPT3(0.5,20,lambda0,phi*pi/180,theta*pi/180,'LINEAL');

lambda_inf = 100;
vol_inf = MAIN_OPT3(0.5,20,lambda_inf,phi*pi/180,theta*pi/180,'LINEAL');


TOL = 1e-2;

lambda_minus = lambda0;
lambda_plus = lambda_inf;
f_minus = vol_0 - Vol;
f_plus = vol_inf - Vol;

while res < TOL
    lambda_new = (lambda_plus - lambda_minus)/2
    vol_new = MAIN_OPT3(0.5,20,lambda_new,phi*pi/180,theta*pi/180,'LINEAL');
    f_new = vol_new - Vol;
       if f_new <= 0
           lambda_minus = lambda_new;
           f_minus = f_new;
       else
           lambda_plus = lambda_new;
           f_plus = f_new;
       end
end
