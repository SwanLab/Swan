function [xNEW, wNEW, nF, ISNEGATIVE,POLYINFO,ISOUT  ] = UpdateCoordinatesPoints(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO)

%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end

ISNEGATIVE = 0 ;
ISOUT = 0 ; 

[PHIk_y,dPHIk_y,POLYINFO ]= EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;

m = length(wNEW);
bNEW = PHIk_y'*wNEW ;
% Equation to be solved
% F =  (b-bNEW)
% Residual
Fk = b-bNEW ;
nF = norm(Fk)  ;
% Computation of the Jacobian matrix
D = [];
for idim = 1:length(dPHIk_y)
    DerPHIt_x = bsxfun(@times,dPHIk_y{idim},wNEW)' ;
    D = [D DerPHIt_x] ;    
end
% Therefore, the Jacobian matrix is finally
D = [D  PHIk_y'] ;
% Updated direction
delta_q = D\Fk ;
%New point
q_kp1 = [xNEW(:);wNEW] +delta_q ;

ndim = length(dPHIk_y); 

wNEW = q_kp1(ndim*m+1:end) ;
%xNEW = [q_kp1(1:m) , q_kp1(m+1:2*m) ] ;

xNEW = [] ; 
INIind= 1 ; 
for idim = 1:ndim    
    FINind = INIind + m -1; 
    xNEW = [xNEW, q_kp1(INIind:FINind) ] ; 
    INIind  = FINind +1 ; 
end


DATA.OnlyCheckIfIsInside = 1;
[ISousideifempty]=     EvaluateBasisFunctionAtX_FEinterp(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ;


if all(wNEW>=0)  && ~isempty(ISousideifempty)
    % disp('All points are admissible !!!!!!!')
  %  disp('')
  
else
    %dbstop('34')
  %  disp('Inadmissible  point ...')
    if any(wNEW<0)   
      %   disp('Some of the weights are negative')
      %   INDLLL = find(wNEW<0); 
      %   wNEG = wNEW(INDLLL) ; 
         
      %   disp(['Negative weights = ',num2str(wNEG'/sum(wNEW)*100),' (%  volume)']) ;       
        ISNEGATIVE = 1 ;
        
    end
    
    if isempty(ISousideifempty)
        ISOUT = 1 ; 
        disp(['Some points are out  of the domain ....'])
    end
    
    
  
end
    



