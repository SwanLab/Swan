function check_derivatives(element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,h_C_0,file_write,iter,x0)

[~,gradient] = equilibrium_update(x0,element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,h_C_0,file_write,iter);
epsilon = 1e-6;
for i = 1:length(x0)

        
        x_plus(:,1) = x0;
        x_minus(:,1) = x0;
        x_plus(i,1) = x_plus(i,1) + epsilon;
        x_minus(i,1) = x_minus(i,1) - epsilon;
        

        [f_plus] = equilibrium_update(x_plus,element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,h_C_0,file_write,iter);
        [f_minus] = equilibrium_update(x_minus,element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,h_C_0,file_write,iter);


        
        dt(i) = (f_plus - f_minus)/(2*epsilon);
        err = norm(dt(i) - gradient(i))/norm(gradient(i))
        error(i) = err;
end

end