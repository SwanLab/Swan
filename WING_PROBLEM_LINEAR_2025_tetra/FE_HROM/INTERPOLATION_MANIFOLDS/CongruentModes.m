function BasisUaligned = CongruentModes(BasisUref,BasisU)
% Inputs:
% REFERENCE MATRIX: BasisUref = [UR1, UR2,...URn]   (where BasisUref'*BasisUref = I)
%   BasisU = [U1,U2...Un]   (where BasisU'*BasisU= I)
% OUTPUT
% BasisUaligned = [V1,V2..Vn], such that
% Vi is a linear combination of U1,U2...Un
% Vi is the most aligned vector of col(BasisU) to the corresponding column
% of BasisUref (Ui)
% JAHO, 11-Dic-2020
%---------------------------------------
if nargin ==0 
    % Example 
    [BasisUref,dummy,dummy2 ]=  SVDT((rand(10,2)));
     BasisU = BasisUref(:,[2,1]) ; 
end 
 

BasisUaligned = [] ;
ORTHOGONAL_compl = BasisU ;
for imodesREF = 1:size(BasisU,2)
    modeREF = BasisUref(:,imodesREF) ;
    % Which is the mode in col(BasisU) most aligned to
    % modeREF  (see hernandez2020multiscale.pdf). Alignment 
    [uu,ss,vv] = SVDT(ORTHOGONAL_compl'*modeREF) ;
    modeALIGNED = ORTHOGONAL_compl*(uu*vv') ;
    BasisUaligned = [BasisUaligned modeALIGNED] ;
    %   [BasisUaligned,SSS,VVV] =  SVDT(BasisUaligned);
    % Orthogonal complement
    ORTHOGONAL_compl = ORTHOGONAL_compl - BasisUaligned*(BasisUaligned'*ORTHOGONAL_compl);
    [ORTHOGONAL_compl,~,~] = SVDT(ORTHOGONAL_compl) ;
end


if nargin == 0
    disp('This matrix should be zero...')
    BasisUaligned-BasisUref
end