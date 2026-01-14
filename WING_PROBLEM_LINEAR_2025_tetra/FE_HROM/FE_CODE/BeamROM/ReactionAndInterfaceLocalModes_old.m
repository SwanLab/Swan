% function [BasisRdef,BasisINT,nBOUNDARY_INTFMODES] = ...
%     ReactionAndInterfaceLocalModes(BasisUdef,BasisRdef,f1,f2,TOL_SINGULAR_VALUES_Hqr,...
%     nBASES_BEAM,DATA_REFMESH,Vrb,M,DISP_BOUND,DATAIN,SinvVal_Udef)
% if nargin == 0
%     load('tmp1.mat')
% end
% ndim = 3;
% f = [f1;f2] ;
% % We have distinguish between "beam" modes and remaining modes (local effects)
% %
% BasisUdef_B = BasisUdef(:,1:nBASES_BEAM.DISPLACEMENTS) ; % Displacements, beam modes
% 
% BasisUdef_L = BasisUdef(:,nBASES_BEAM.DISPLACEMENTS+1:end) ; % Displacements, local effects
% SingVal_Udef_L = SinvVal_Udef(nBASES_BEAM.DISPLACEMENTS+1:end)  ; % Associated singular values
% BasisRdef_L = BasisRdef(:,nBASES_BEAM.REACTIONS+1:end) ;     % Reactions, local effects  
% % What to do with these modes  ?  
% %  Reference matrix 
%  
% % Candidates to be interface displacement modes
% % ----------------------------------------------
%  [nBOUNDARY_INTFMODES,BasisREFERENCE] = CandidatesInterfaceModes(BasisUdef_L,f1,f2,DISP_BOUND,SingVal_Udef_L,...
%      DATAIN,BasisUdef_B) ; 
% 
% 
%  
% % BasisUdef_Lf1 = BasisUdef_L(f_reference,:) ;
% % BasisRdef_Lf1 = BasisRdef_L(f_reference,:) ;
% 
% BasisUdef_Lf = BasisUdef_L(f,:) ;
% BasisRdef_Lf = BasisRdef_L(f,:) ;
% 
% % Determining number of reaction modes
% %b -----------------------------------
% imode = 1;
% 
% MODES_INCLUDE =[] ;
% IMODES_INCLUDE = [] ;
% while  imode <=size(BasisRdef_Lf,2)
%     NEW_MODES = [MODES_INCLUDE,BasisRdef_Lf(:,imode) ] ;
%     HqrT = NEW_MODES'*BasisUdef_Lf;
%     SSVAL = svd(HqrT) ;
%     if imode == 1 
%         ratioSV =1 ; % SSVAL(1) ; 
%     else
%     ratioSV = SSVAL(end)/SSVAL(end-1) ;
%     end
%     if ratioSV >= TOL_SINGULAR_VALUES_Hqr
%         MODES_INCLUDE = NEW_MODES ;
%         IMODES_INCLUDE(end+1) = imode ;
%     end
%     imode = imode + 1;
%     
%     
% end
% 
% nmodesR = imode-1;
% BasisRdef_L = BasisRdef_L(:,IMODES_INCLUDE);
% BasisRdef = [BasisRdef(:,1:nBASES_BEAM.REACTIONS) BasisRdef_L] ;
% 
% DATAOUT.BasisRdef = BasisRdef ;
% % Option 1
% % --------
% 
% 
% 
% %USE_MASS_MATRIX = 1 ; 
% %if USE_MASS_MATRIX == 0 
% %    M = speye(size(M)) ; 
% %end 
% 
% PG = (Vrb'*M*Vrb)  ;
% 
% 
% 
% %     BasisINTdef = [] ;
% %
% %     for imode = 1:size(BasisRdef_L,2)
% %
% %     BasisINTdef = BasisRdef_L(f1,:)  - Vrb*(PG\(Vrb'*M*BasisRdef_L(f1,:))) ;
% %
% %     end
% 
% %WHICH_MODES_USE = 1 ; 
% %if WHICH_MODES_USE == 0
% %    BasisREFERENCE = BasisRdef_L(f1,:) ; 
% %elseif WHICH_MODES_USE == 1
% %end
% 
% MODES_INCLUDE =[] ;
% IMODES_INCLUDE = [] ;
% imode = 1;
% 
% % The number of interface modes cannot be greater than the number of
% % reaction modes
% while  imode <=size(BasisRdef_L,2)
% %    if USE_REACTIONS == 1
%     newINTFmode = BasisREFERENCE(:,imode)  - Vrb*(PG\(Vrb'*M*BasisREFERENCE(:,imode))) ;
%  %   else
%    %     newINTFmode = BasisRdef_L(f1,imode)  - Vrb*(PG\(Vrb'*M*BasisRdef_L(f1,imode))) ;
%   %  end
%     NEW_MODES = [MODES_INCLUDE,newINTFmode ] ;
%     COV_f1 = NEW_MODES'*BasisRdef_L(f1,:);
%     COV_f2 = NEW_MODES'*BasisRdef_L(f2,:);
%     SSVAL_f1 = svd(COV_f1) ;
%     SSVAL_f2 = svd(COV_f2) ;
%     if imode == 1
%       %  if SSVAL_f1(1) > 1e-10 &  
%         ratioSV_f1 = 1  ; 
%         ratioSV_f2 = 1  ; 
%        % end
%     else
%     ratioSV_f1 = SSVAL_f1(end)/SSVAL_f1(end-1) ;
%     ratioSV_f2 = SSVAL_f2(end)/SSVAL_f2(end-1) ;
%     end
%     if ratioSV_f1 >= TOL_SINGULAR_VALUES_Hqr && ratioSV_f2 >= TOL_SINGULAR_VALUES_Hqr
%         MODES_INCLUDE = NEW_MODES ;
%         IMODES_INCLUDE(end+1) = imode ;
%     end
%     imode = imode + 1;
%     
%     
% end
% 
% 
% BasisINTdef = MODES_INCLUDE ;
% 
% BasisINT = [Vrb,BasisINTdef] ;