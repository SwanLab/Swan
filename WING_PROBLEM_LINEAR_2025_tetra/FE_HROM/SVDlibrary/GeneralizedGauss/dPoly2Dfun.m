function [Xf_x Xf_y ]= dPoly2Dfun(P,x,y) ;

if nargin ==0
    P = 1
    x = sqrt(0.5)
    y = sqrt(0.5)
end



if mod(P,2)==0
    nmo = (P+1)*(P+2)/2 ;
else
    nmo = (P)*(P+1)/2 + (P+1) ;
end

Xf_x= zeros(length(x),nmo) ;
iacum = 1 ;
for i = 0:P
    for j=0:P
        if i+j<=P
            if i >0
                Xf_x(:,iacum) = i*(x.^(i-1)).*(y.^j) ;
             
            end
               iacum = iacum + 1;
        end
    end
end

Xf_y= zeros(length(x),nmo) ;
iacum = 1 ;
for i = 0:P
    for j=0:P
        if i+j<=P
            if j>0
                Xf_y(:,iacum) = j*(x.^i).*(y.^(j-1)) ;
              
            end
              iacum = iacum + 1;
        end
    end
end

end