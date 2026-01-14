function [BasisUnew,SingularValuesNewDisp,TEXT,BasisRnew,SingularValuesNewREAC ]= ReactionDispModesEqual...
    (BasisUdef_L,BasisRdef_L,f,SingVal_Udef,TEXT,DATAIN)
% See LATEX  ReactionDispModesEqual.tex
if nargin == 0
    load('tmp1.mat')
elseif nargin == 5
    DATAIN = [] ;
end


if size(BasisUdef_L,2) == size(BasisRdef_L,2)
    BasisUnew = BasisUdef_L ;
    % TEXT = {} ;
    SingularValuesNewDisp = SingVal_Udef ;
    TEXT{end+1} = ['Number react modes = Number displ. modes = ',num2str(size(BasisUnew,2))] ;
    U = BasisUdef_L ; 
    F = BasisRdef_L ; 
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%    warning('Remove that...')
 %   SingVal_Udef = ones(size(SingVal_Udef)) ;  % remove
    Y = bsxfun(@times,BasisUdef_L',SingVal_Udef)' ;  % Snaphots deformational displacements
    z = BasisRdef_L(f,:)  ; % Basis matrix for forces
    [y,S,V ]= SVDT(Y(f,:),0) ; % Basis matrix boundary deformational disp.
    SV = bsxfun(@times,V',S)' ;
    % Matrix A (least squares)
    % --------
    A = y'*z ;
    [UA,SA,VA] = SVDT(A,0) ;
    A = UA*VA';
    % Coefficients
    SVA = SV*A ;
    % New matrix
    BasisUnew =  BasisUdef_L*SVA ;
    [BasisUnew,SingularValuesNewDisp ]=SVDT(BasisUnew,0) ;
    
    
    %%% COMPUTE INTERSECTION --- Curiosity 2-OCT-2019
    % --------------------------------------------------------
    % [yA,~,~] = SVDT(y,0) ;
    % [yB,~,~] = SVDT(z,0) ;
    % [ZZ,SS,~] = SVDT(yA'*yB,0) ;
    % aaaa =acos(SS)*180/pi ;
    
    %%%%%
    
    if length(SingularValuesNewDisp) ~=size(z,2)
        
        disp(['Check this routine. Number of displacement modes less than number of reactions'])
        disp(['Set a lower tolerance for the SVD of def. displacements'])
        error(' ')
    end
    
    U = BasisUnew(f,:) ;
    F= BasisRdef_L(f,:) ;
    
    H = U'*F;
    TEXT = {} ;
    SSS = svd(H) ;
    TEXT{end+1} = ['Number react modes = Number displ. modes = ',num2str(size(z,2))] ;
    TEXT{end+1} = ['svd(H,end)/svd(H,1) ',num2str(SSS(end)/SSS(1))] ;
    TEXT{end+1} = ['svd(H,end)/svd(H,end-1) ',num2str(SSS(end)/SSS(end-1))] ;
    
    
end


% MODIFICATION 13-MAY-2020.Issues encountered  in
%  /home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/NEW_IDEAS/DOCS/RAUL_BRAVO_THESHIS/TowardsMinimumNPoints/PLASbeam/Consistency_nonl.m
% Apparently, when the angles formed by the reaction and displacements are
% close to 90 degrees, the matrix H becomes ill-posed. In what follows we
% propose a procedure that attemtps to alleviate this issue

DATAIN = DefaultField(DATAIN,'FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES',0) ;
DATAIN = DefaultField(DATAIN,'ANGLE_FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES',89) ; % Degrees
TOL_SVD = 0 ;
    ANGLE = DATAIN.ANGLE_FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES ;
    [RA,RB,ANGLES ]= AlignedSubspaces(BasisUnew,BasisRdef_L,TOL_SVD,ANGLE) ;
    
    TEXT{end+1} = ['Checking angles between subspace REACT and DEF = ',num2str(ANGLES')] ; 
    
if DATAIN.FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES ==1
    
    BasisUnew = RA;
    BasisRnew = RB ;
    SingularValuesNewDisp = ones(size(RA,2),1) ;
    SingularValuesNewREAC = ones(size(RB,2),1) ;
else
    BasisRnew = BasisRdef_L ;
    SingularValuesNewREAC = [] ;
   
end


