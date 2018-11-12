function stop = outfun(x,optimValues,state)
stop = false;
global history 
   switch state
       case 'init'

       case 'iter'
           
          
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
           history.constr = [history.constr; optimValues.constrviolation];
           
           figure(1)
           plot([1:optimValues.iteration+1],history.fval)


       case 'done'
           hold off
       otherwise
   end
end