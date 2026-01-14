function [nREDcoor,Autoencoder] = ConstructAutoencoder(PhiLIN,PhiNON,DATA_interp,qINF,VVt,UU,SS,qSUP)

if nargin == 0
    load('tmp.mat')
end

BasisU = [PhiLIN,PhiNON] ;
if ~isempty(qSUP)     
    DATA_interp = DefaultField(DATA_interp,'SubSampling_qINELASTIC_master',false) ;    
    if ~DATA_interp.SubSampling_qINELASTIC_master        
        switch  DATA_interp.METHOD_SELECT_REFERENCE_MODE
            case  {'MAXIMUM_DISPLACEMENT_LOCATION','DISPLACEMENT_ONE_NODE'}
        [DATA_evaluateTAU_and_DER,nREDcoor]= BsplinesLeastSquares_MAXLOC(DATA_interp, qINF, VVt, UU, SS) ; 
%             case 'DISPLACEMENT_ONE_NODE'
%                  [DATA_evaluateTAU_and_DER,nREDcoor]= BsplinesLeastSquares_DISP1node(DATA_interp, qINF, VVt, UU, SS) ; 
            otherwise
            [DATA_evaluateTAU_and_DER,nREDcoor]= BsplinesLeastSquares_fast(DATA_interp, qINF, VVt, UU, SS) ;     
        end
        
    else
        [DATA_evaluateTAU_and_DER,DATA_interp,nREDcoor]= BsplinesSUBSAMPL_large(qINF,qSUP,DATA_interp ) ;
    end    
else    %[DATA_evaluateTAU_and_DER,nREDcoor]= StandardROM(DATA_interp,BasisU ) ;    
    nREDcoor = size(BasisU,2); 
    DATA_evaluateTAU_and_DER.nameFunctionEvaluate = 'tauFUN_identity';    
end