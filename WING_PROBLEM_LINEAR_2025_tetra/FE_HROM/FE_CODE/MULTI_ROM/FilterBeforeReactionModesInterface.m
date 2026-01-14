function [BasisRdef1,BasisRdef2] = FilterBeforeReactionModesInterface(DATAIN,BasisRdef,f1,f2)

if ~isempty(DATAIN.FILTER_OUT_REACTION_MODES_RATIO_SINGULAR_VALUES)
    [UR1,SR1,VR1] = SVDT(BasisRdef(f1,:)) ;
    rS1 = SR1(2:end)./SR1(1:end-1) ;
    nS1 = (find(rS1 < DATAIN.FILTER_OUT_REACTION_MODES_RATIO_SINGULAR_VALUES)) ;
    if isempty(nS1)
        nS1 = size(BasisRdef,2) ;
    end
    [UR2,SR2,VR2] = SVDT(BasisRdef(f2,:)) ;
    rS2 = SR2(2:end)./SR2(1:end-1) ;
    nS2 = (find(rS2 < DATAIN.FILTER_OUT_REACTION_MODES_RATIO_SINGULAR_VALUES)) ;
    if isempty(nS2)
        nS2 = size(BasisRdef,2) ;
        
        
    end
    
    nS = min(nS1,nS2) ;
    BasisRdef1  =UR1(:,1:nS) ;
    BasisRdef2  =UR2(:,1:nS) ;
else
    BasisRdef2 = BasisRdef(f2,:) ;
    BasisRdef1 = BasisRdef(f1,:) ;
    
end