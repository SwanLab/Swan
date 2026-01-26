function [BasisINT,RotationMatrixLOC,TEXTP] =  ReactionAndInterfaceLocalModes_RVE(BasisUdef,BasisRdef,f1,f2,...
    TOL_SINGULAR_VALUES_Hqr,...
    Vrb,M,DATAIN,SinvVal_Udef,SinvVal_Rdef,iface,TEXTP)
% Copy of ReactionAndInterfaceLocalModes_new  (this one was made for beam structures)
%  JAHO
if nargin == 0
    load('tmp.mat')
end

f = [f1;f2] ;
RotationMatrixLOC = [] ;   

% Candidates to be interface fluctuation displacement modes. This function
% also orthogonalize Vrb with respect to M
% ----------------------------------------------
[VrbORTH,VdefORTH] = CandidatesInterfaceModesRVE(BasisUdef,f1,f2,M,Vrb,DATAIN) ;
nmodesRB = size(VrbORTH,2) ;

MODES_INCLUDE =[] ;
IMODES_INCLUDE = [] ;
imode = 1;


DATAIN = DefaultField(DATAIN,'DeformationalInterfaceModes_AlignmentMethod',0) ;
DATAIN = DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-3) ;

if DATAIN.DeformationalInterfaceModes_AlignmentMethod == 1
     % warning('This implementation is not reliable')
    
   [ BasisINT,RotationMatrixLOC,TEXTP ]=  DeformModesInterface_AlignmentMethod(BasisRdef,f1,f2,...
        Vrb,M,DATAIN,BasisUdef,SinvVal_Udef,SinvVal_Rdef,iface,TEXTP) ;
    % Rotated Reactions ---> Reactions expressed in the ref. system
    % attached to the boundary interfaces, 9-Apr-2019
    
else
    Vcandidate = [VrbORTH,VdefORTH] ;
    
    
    % The number of interface modes cannot be greater than the number of
    % reaction modes
    while  imode <=size(BasisRdef,2) &&  imode <=size(Vcandidate,2)
        %    if USE_REACTIONS == 1
        newINTFmode = Vcandidate(:,imode) ;
        %   else
        %     newINTFmode = BasisRdef(f1,imode)  - Vrb*(PG\(Vrb'*M*BasisRdef(f1,imode))) ;
        %  end
        NEW_MODES = [MODES_INCLUDE,newINTFmode ] ;
        COV_f1 = NEW_MODES'*BasisRdef(f1,:);
        COV_f2 = NEW_MODES'*BasisRdef(f2,:);
        SSVAL_f1 = svd(COV_f1) ;
        SSVAL_f2 = svd(COV_f2) ;
        if imode == 1
            %  if SSVAL_f1(1) > 1e-10 &
            ratioSV_f1 = 1  ;
            ratioSV_f2 = 1  ;
            % end
        else
            ratioSV_f1 = SSVAL_f1(end)/SSVAL_f1(end-1) ;
            ratioSV_f2 = SSVAL_f2(end)/SSVAL_f2(end-1) ;
        end
        if ratioSV_f1 >= TOL_SINGULAR_VALUES_Hqr && ratioSV_f2 >= TOL_SINGULAR_VALUES_Hqr
            MODES_INCLUDE = NEW_MODES ;
            IMODES_INCLUDE(end+1) = imode ;
        else
            if nmodesRB> imode
                error('Enrich the training set. Reaction modes do not form a consistent set.')
            end
        end
        imode = imode + 1;
        
        
    end
    
%normVrb  =norm(Vrb,'fro') ; 
%normDEF = norm(MODES_INCLUDE(:,nmodesRB+1:end),'fro')  ; 
%Vdef = MODES_INCLUDE(:,nmodesRB+1:end)*normVrb/normDEF ; 
  %  BasisINT = [Vrb,Vdef] ;

    BasisINT = [Vrb, MODES_INCLUDE(:,nmodesRB+1:end)] ;
   % BasisINT = [VrbORTH, MODES_INCLUDE(:,nmodesRB+1:end)] ;

% 
% % Verification
% [UVer1,Sver1,Vver1 ] = SVDT(BasisINT'*BasisRdef(f1,:)) ;
% [UVer2,Sver2,Vver2 ] = SVDT(BasisINT'*BasisRdef(f2,:)) ;
    
end


