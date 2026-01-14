function E = ConvertFgrad_greenlg(GRADuMACRO) 

if size(GRADuMACRO,1) == 3
    
%     Egreenl_T = zeros(3,3) ; 
%     Egreenl_T(1,1) = Egreenl(1) ; 
%     Egreenl_T(2,2) = Egreenl(2) ; 
%     Egreenl_T(3,3) = Egreenl(3) ; 
%     Egreenl_T(2,3) = Egreenl(4) ; 
%     Egreenl_T(1,3) = Egreenl(5) ; 
%     Egreenl_T(1,2) = Egreenl(6) ; 
%     Egreenl_T(3,2) = Egreenl(4) ; 
%     Egreenl_T(3,1) = Egreenl(5) ; 
%     Egreenl_T(2,1) = Egreenl(6) ; 
%     
%     FtF = 2*Egreenl_T + eye(3) ;  
%     
%     % Cholesky decomposition 
%     Ft = chol(FtF) ;
%     F = Ft' ; 
%     [U,S,V] =svd(F) ; 
%     
%     % Polar decomposition 
%     Fdef = V*S*V' ; 
    
    %Fgrad  =F;  
    Fgrad = GRADuMACRO + eye(3) ; 
    Egreenl_T = 0.5*(Fgrad'*Fgrad - eye(3)) ; 
    
    EgreenEFF = zeros(6,1) ; 
    
      EgreenEFF(1) = Egreenl_T(1,1) ;
     EgreenEFF(2) = Egreenl_T(2,2) ;
     EgreenEFF(3) = Egreenl_T(3,3)  ;
     EgreenEFF(4) = Egreenl_T(2,3) ;
     EgreenEFF(5) = Egreenl_T(1,3)  ;
     EgreenEFF(6) =  Egreenl_T(1,2)  ;
   
    
    disp(['-------------------------------------------------'])
    disp(['Effective Green-Lagrange strain tensor  '])
    disp(['E = ',num2str(EgreenEFF')])
     disp(['-------------------------------------------------'])
    
else
    error('Option not implemented yet')
end 
