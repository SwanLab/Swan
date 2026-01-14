function s = Voig2MatSTRAIN(s) ;


if length(s) == 6

s = [s(1)  0.5*s(6)  0.5*s(5)
     0.5*s(6)  s(2)  0.5*s(4)
     0.5*s(5)  0.5*s(4)  s(3)] ; 
else
    
    s = [s(1)  0.5*s(3)  
     0.5*s(3)  s(2) ] ;
    
end

