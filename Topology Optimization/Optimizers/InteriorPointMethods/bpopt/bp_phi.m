function [ph] = bp_phi(bp,x,xL,xU,s,bL,bU)
    ph = bp_obj(bp,x);
    n = size(x,2);
    m = size(bL,2);
    for i = 1:n,
       ph = ph - bp.mu * (log(x(i)-xL(i)) + log(xU(i)-x(i)));
    end
    j = 0;
    for i = 1:m,
       if(bU(i)>bL(i)),
          j = j + 1;
          ph = ph - bp.mu * (log(s(j)-bL(i)) + log(bU(i)-s(j)));
       end
    end
