function [xNEW,wNEW,ISNEGATIVE,ISOUT] =  UpdateSolutionWP_FORCED(xNEW,wNEW,delta_q,DATA,VAR_SMOOTH_FE,POLYINFO)

if nargin == 0
    load('tmp.mat')
end

ISNEGATIVE = 0 ;
ISOUT = 0 ; 


%New point
xOLD = xNEW ; 
q_kp1 = [xNEW(:);wNEW] +delta_q ;

ndim = size(xNEW,2) ;  
m = length(wNEW);
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
[ISousideifempty,~,POLYINFO_NEW]=     EvaluateBasisFunctionAtX_LARGE(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ;

if  ~isempty(POLYINFO_NEW.ListPointsOutside)
    disp([num2str(length(POLYINFO_NEW.ListPointsOutside)), ' points are being returned to the domain ',])
    
    xNEW(POLYINFO_NEW.ListPointsOutside,:) = xOLD(POLYINFO_NEW.ListPointsOutside,:)  ; 
    ISousideifempty = 0 ; 
end


ElementsBefore = POLYINFO.ELEMENTS_CONTAINING_xNEW; 
ElementsNew = POLYINFO_NEW.ELEMENTS_CONTAINING_xNEW ; 
DiffElements = ElementsNew-ElementsBefore ; 
AAA =find(DiffElements ~=0) ;
if isempty(AAA)
else
  %  disp(['Interelement Jump']) 
    ElementsJump_before_after = [ElementsBefore(AAA),ElementsNew(AAA)] ;
    
    ElementsJump_before_after_all = unique(ElementsJump_before_after(:)) ; 
    ElementsJump_before_after_all = ElementsJump_before_after_all(ElementsJump_before_after_all>0); 
    
end


NNEGATIVES  = length(find(wNEW<0)) ; 


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
         disp(['Found ',num2str(length(find(wNEW<0))),' negative weights']); 
        if NNEGATIVES > DATA.THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS
            ISOUT = 1; 
            disp(['Number of negative points above the imposed threshold'])
        end
        
    end
    
    if isempty(ISousideifempty)
        ISOUT = 1 ; 
        disp(['Some points are out  of the domain ....'])
    end
    
    
  
end
    



