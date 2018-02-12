function check_derivatives_constr(dim,element,problembsc,coordinates,x0)

[~,gradient] = call_constraint(x0,dim,element,problembsc,coordinates);
epsilon = 1e-6;
for i = 1:length(x0)

        
        x_plus(:,1) = x0;
        x_minus(:,1) = x0;
        x_plus(i,1) = x_plus(i,1) + epsilon;
        x_minus(i,1) = x_minus(i,1) - epsilon;
        

        [f_plus] = call_constraint(x_plus,dim,element,problembsc,coordinates);
        [f_minus] = call_constraint(x_minus,dim,element,problembsc,coordinates);


        
        dt(i) = (f_plus - f_minus)/(2*epsilon);
        err = norm(dt(i) - gradient(i))/norm(gradient(i))
        error(i) = err;
end

end