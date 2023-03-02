function [ac] = bp_accept(bp,x,xa,xL,xU,filter)
% parameters
filter_on = false;
better_on = false;
n = size(x,2);
gamma_th = 10^-5;
gamma_phi = 10^-5;

% check for violations of variable constraints
% shouldn't need to check these
viol = 0;
for i = 1:n,
    if(xa(i)<xL(i)),
        viol = viol + 1;
    end
    if(xa(i)>xU(i)),
        viol = viol + 1;
    end
end

if (viol==0),
  % either condition better compared to last iteration
  if (better_on),
     better = (bp_phi(bp,xa,xL,xU) <= bp_phi(bp,x,xL,xU) - gamma_phi*bp_theta(bp,x)) || (bp_theta(bp,xa) <= (1-gamma_th)* bp_theta(bp,x));
  else
     better = true;
  end
     
  % apply filter to determine whether to accept
  if (better),
      ac = true;
      if (filter_on),
         nf = size(filter,1);
         for i = 1:nf,
             if (((1-gamma_th)*bp_theta(bp,xa) > filter(i,1)) || (bp_phi(bp,xa,xL,xU)-gamma_phi*bp_theta(bp,xa) > filter(i,2))),
                 ac = false;
             end
         end
      end
  else
      ac = false;
  end
else
    ac = false;
end
