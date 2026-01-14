function [xNEW,wNEW,ISNEGATIVE,ISOUT] =  UpdateSolutionWP_NOTFORCED(xNEW,wNEW,delta_q,DATA,VAR_SMOOTH_FE,POLYINFO)

if nargin == 0
    load('tmp1.mat')
end

ISNEGATIVE = 0 ;
ISOUT = 0 ; 


%New point
ndim = size(xNEW,2) ;  
m = length(wNEW);

if DATA.LINE_SEARCH_AVOID_NEGATIVE_WEIGTHS ==0
    STEP_SIZE = 1; 
else
    error('This option turned out to be unreliable')
    DELTA_qw = delta_q(ndim*m+1:end) ;
    wNEW_tentative = DELTA_qw+wNEW ; 
    IND_NEGATIVES = find(wNEW_tentative <0 ) ; 
  %  [minW,indMIN] = min(wNEW_tentative) ; 
    if isempty(IND_NEGATIVES)  
        STEP_SIZE = 1; 
    else
       
        delta_critic_w = DELTA_qw(IND_NEGATIVES) ;   % These are  the entries of the incremental weight
        % that render some weights negative (dwn)       
        delta_for_zero = -wNEW(IND_NEGATIVES) ;  % These the value of the increment for the weights 
        % to remain positive  
        STEP_SIZE_all = delta_for_zero./delta_critic_w ; 
        STEP_SIZE = min(STEP_SIZE_all); 
        
        if  STEP_SIZE == 0
            % Ignoring this 
            STEP_SIZE_w = ones(size(DELTA_qw)); 
            STEP_SIZE_w(IND_NEGATIVES) =   STEP_SIZE_all ; 
            STEP_SIZE = ones(size(delta_q)) ; 
            STEP_SIZE(ndim*m+1:end) = STEP_SIZE_w ; 
        end
        
        
    end
end


q_kp1 = [xNEW(:);wNEW] +STEP_SIZE.*delta_q ;


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
ElementsBefore = POLYINFO.ELEMENTS_CONTAINING_xNEW; 
ElementsNew = POLYINFO_NEW.ELEMENTS_CONTAINING_xNEW ; 
DiffElements = ElementsNew-ElementsBefore ; 
AAA =find(DiffElements ~=0) ;
if isempty(AAA)
else
   % disp(['Interelement Jump']) 
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
    



