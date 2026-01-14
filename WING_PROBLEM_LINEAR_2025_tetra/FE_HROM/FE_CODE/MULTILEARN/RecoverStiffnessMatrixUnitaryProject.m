function  K = RecoverStiffnessMatrixUnitaryProject(DATAONLINE)

if nargin == 0
    load('tmp.mat')
end


 PROJECT_LOADS = DATAONLINE.PROJECT_LOADS{1} ; 
 eval(PROJECT_LOADS) ; 
DATA.INPUTDATAfile = PROJECT_LOADS  ;
DATA.NOCALCULATE_DISPLACEMENTS = 1 ;
% Calling Finite Element elastostatic program (but only for computing external forces)
DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
load(DATAOUT.nameWORKSPACE,'K') ;

% [IDX D]= knnsearch(COOR,COORref) ;
% 
% if any(abs(D) >1e-16 )
%     dbstop('51')
%     error('Non-conforming meshes')
% end
% IDXdofs = Nod2DOF(IDX,size(COOR,2)) ;
% % All properties are set in terms of the numbering of reference mesh
% %  dbstop('55')
%      load(DATAOUT.nameWORKSPACE,'Cglo') ;
%     KdomRED = {} ;
%     for itype = 1:length(BdomRED)
%         nstrain = size(BdomRED{1},1)/length(Wdom) ;
%         if  DATAINM.CUBATURE.ACTIVE == 1
%             setIndices =  small2large(setPoints{itype},nstrain) ;
%             % Cglo is multiplied by Wdom(setPoints). Accordingly, we
%             % define
%             rW = WdomRED{itype}./Wdom(setPoints{itype}) ;
%             % And make
%             rwDIAG = CompWeightDiag(rW,nstrain)  ;
%             KdomRED{itype} = (rwDIAG*BdomRED{itype}(setIndices,:))'*(Cglo(setIndices,setIndices)*BdomRED{itype}(setIndices,:)) ;
%             
%         else
%             KdomRED{itype} =   BdomRED{itype}'*(Cglo*BdomRED{itype}) ;
%         end
%     end
%     
%  
%  