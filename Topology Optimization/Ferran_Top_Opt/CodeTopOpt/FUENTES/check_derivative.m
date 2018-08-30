function error = check_derivative(x0,fun)

[f0,gradient] = fun(x0);
epsilon = 1e-6;
for i = 1:length(x0)

        
        x_plus(:,1) = x0;
        x_minus(:,1) = x0;
        x_plus(i,1) = x_plus(i,1) ;%+ epsilon;
        x_minus(i,1) = x_minus(i,1) - epsilon;
        

        [f_plus] = fun(x_plus);
        [f_minus] = fun(x_minus);

        
        dt(i) = (f_plus - f_minus)/epsilon;%(2*epsilon);
        err = norm(dt(i) - gradient(i))/norm(dt(i));
        error(i,1) = err
end

end