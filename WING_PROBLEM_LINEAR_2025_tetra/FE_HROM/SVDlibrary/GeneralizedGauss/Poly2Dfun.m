function Xf = Poly2Dfun(P,x,y) ;

 
 
%dbstop('6')
if mod(P,2)==0
    nmo = (P+1)*(P+2)/2 ; 
else
     nmo = (P)*(P+1)/2 + (P+1) ; 
end

Xf = zeros(length(x),nmo) ; 
iacum = 1 ; 
for i = 0:P
    for j=0:P
        if i+j<=P
        Xf(:,iacum) = (x.^i).*(y.^j) ;
        iacum = iacum + 1;
        end
    end
end