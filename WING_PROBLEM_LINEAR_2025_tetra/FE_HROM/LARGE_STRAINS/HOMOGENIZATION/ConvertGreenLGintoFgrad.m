function [Fgrad,GradU ]= ConvertGreenLGintoFgrad(Egreenl)

if nargin == 0
    load('tmp.mat')
end

if length(Egreenl) == 6
    
    Egreenl_T = zeros(3,3) ;
    Egreenl_T(1,1) = Egreenl(1) ;
    Egreenl_T(2,2) = Egreenl(2) ;
    Egreenl_T(3,3) = Egreenl(3) ;
    Egreenl_T(2,3) = 0.5*Egreenl(4) ;
    Egreenl_T(1,3) = 0.5*Egreenl(5) ;
    Egreenl_T(1,2) = 0.5*Egreenl(6) ;
    Egreenl_T(3,2) = 0.5*Egreenl(4) ;
    Egreenl_T(3,1) = 0.5*Egreenl(5) ;
    Egreenl_T(2,1) = 0.5*Egreenl(6) ;
    
    FtF = 2*Egreenl_T + eye(3) ;
    
    % Cholesky decomposition
    F = chol(FtF) ;
    %S  F = Ft' ;
    [U,S,V] =svd(F) ;
    
    % Polar decomposition
    Fdef = V*S*V' ;
    
    Fgrad  =F;
    GradU = F-eye(3) ;
    
    Egreenl_T = 0.5*(Fgrad'*Fgrad - eye(3)) ;
    
    EgreenEFF = zeros(6,1) ;
    
    EgreenEFF(1) = Egreenl_T(1,1) ;
    EgreenEFF(2) = Egreenl_T(2,2) ;
    EgreenEFF(3) = Egreenl_T(3,3)  ;
    EgreenEFF(4) = Egreenl_T(2,3) ;
    EgreenEFF(5) = Egreenl_T(1,3)  ;
    EgreenEFF(6) =  Egreenl_T(1,2)  ;
    
    
    disp(['-------------------------------------------------']);
    disp(['Effective Green-Lagrange strain tensor  ']);
    disp(['E = ',num2str(EgreenEFF')]);
    disp(['-------------------------------------------------']);
    
else
    % This is only for tHE PLANE STRAIN CASE 
    Egreenl_T = zeros(2,2) ;
    Egreenl_T(1,1) = Egreenl(1) ;
    Egreenl_T(2,2) = Egreenl(2) ;
      Egreenl_T(1,2) =Egreenl(3);
    Egreenl_T(2,1) = Egreenl(3)  ;
    
    FtF = 2*Egreenl_T + eye(2) ;
    
    % Cholesky decomposition
   % F = chol(FtF) ;
     F = cholcov(FtF) ; 
     if isempty(F)
         error('It is  not possible to determine F')
     end
     if  size(F,1) ~= size(FtF,1)
        Fnew= zeros(size(FtF)) ;         
        nrows = size(Fnew,1)-size(F,1)+1: size(Fnew,1) ;         
        Fnew(nrows,:) = F ;         
        F = Fnew  ;          
     end
     
     
    %S  F = Ft' ;
   % [U,S,V] =svd(F) ;
    
    % Polar decomposition
    %Fdef = S*V' ;
    
    Fgrad  = F ;
    GradU = F-eye(2) ;
    
    
    
    Egreenl_T = 0.5*(Fgrad'*Fgrad - eye(2)) ;
    
    EgreenEFF = zeros(3,1) ;
    
    EgreenEFF(1) = Egreenl_T(1,1) ;
    EgreenEFF(2) = Egreenl_T(2,2) ;
    EgreenEFF(3) = Egreenl_T(1,2)  ;
%     EgreenEFF(4) = Egreenl_T(2,3) ;
%     EgreenEFF(5) = Egreenl_T(1,3)  ;
%     EgreenEFF(6) =  Egreenl_T(1,2)  ;
    
    
    disp(['-------------------------------------------------']);
    disp(['Effective Green-Lagrange strain tensor  ']);
    disp(['E = ',num2str(EgreenEFF')]);
    disp(['-------------------------------------------------']);
    
end
