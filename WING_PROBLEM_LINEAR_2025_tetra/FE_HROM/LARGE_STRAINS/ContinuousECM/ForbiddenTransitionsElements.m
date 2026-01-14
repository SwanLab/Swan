function [VARCnew,ListForbiddenTransitions]  = ForbiddenTransitionsElements(POLYINFO,DATA,POLYINFO_NEW,VARCnew)

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx

ElementsBefore = POLYINFO.ELEMENTS_CONTAINING_xNEW;
ElementsNew = POLYINFO_NEW.ELEMENTS_CONTAINING_xNEW ;
DiffElements = ElementsNew-ElementsBefore ;
AAA =find(DiffElements ~=0) ;
ListForbiddenTransitions = [] ;
if isempty(AAA)
    VARCnew.ListElementsInTransitionINNERloop{end+1} = {} ;
    
else
    
    ElementsJump_before_after = [ElementsBefore(AAA),ElementsNew(AAA)] ;
    
    ElementsJump_before_after_all = unique(ElementsJump_before_after(:)) ;
    ElementsJump_before_after_all = ElementsJump_before_after_all(ElementsJump_before_after_all>0);
    
    if  length(ElementsJump_before_after_all)>1
     %   disp(['Interelement Jump: ',num2str(ElementsJump_before_after_all')])
        INZZ = find(ElementsNew(AAA)) ;
        ElementsInTransition = unique(ElementsJump_before_after(INZZ,:)) ;
        VARCnew.ListElementsInTransitionINNERloop{end+1} = ElementsInTransition(:) ;
        if ~isempty(DATA.EXCLUDE_ELEMENT_TRANSITION_LIST) && ~ischar(DATA.EXCLUDE_ELEMENT_TRANSITION_LIST)
            for itransition = 1:size(ElementsJump_before_after,1)
                isBefore =  ismember(ElementsJump_before_after(itransition,1),DATA.EXCLUDE_ELEMENT_TRANSITION_LIST) ;
                isAfter=  ismember(ElementsJump_before_after(itransition,2),DATA.EXCLUDE_ELEMENT_TRANSITION_LIST) ;
                if isBefore && isAfter
                    ipoint_to_freeze = AAA(itransition) ;
                    disp(['Point =',num2str(ipoint_to_freeze), ' is not allowed to move from element ',...
                        num2str(ElementsJump_before_after(itransition,1)),' to element ',num2str(ElementsJump_before_after(itransition,2)),...
                        '(freezing its position ...)'])
                    ListForbiddenTransitions = [ListForbiddenTransitions;ipoint_to_freeze] ;
                end
            end
        end
    else
        VARCnew.ListElementsInTransitionINNERloop{end+1} = {} ;
    end
end