function [BasisUnew]= ReactionDispModesEqual2022(AsnapDISPdef,PsiDEFf,f,DATAROM)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DOMAINdecompositionGEN/TruncationAnalysis.mlx
if nargin == 0
    load('tmp3.mat')
elseif nargin == 3
    DATAROM = [] ;
end

[BasisUdef_L,SingVal_Udef,~] = SVDT(AsnapDISPdef) ;

%INCLUDE_SINGULAR_VALUES = 2 ;
% See discussion in  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DOMAINdecompositionGEN/TruncationAnalysis.mlx
%if INCLUDE_SINGULAR_VALUES == 1
%    Y = bsxfun(@times,BasisUdef_L',SingVal_Udef)' ;  % Snaphots deformational displacements
%elseif INCLUDE_SINGULAR_VALUES == 0
%   Y = BasisUdef_L ;
%elseif INCLUDE_SINGULAR_VALUES == 2
Y = AsnapDISPdef ;
%end
%[z,~,~ ]=SVDT(PsiDEFf) ; % Basis matrix for forces

z = PsiDEFf ; 

[y,S,V ]= SVDT(Y(f,:),0) ; % Basis matrix boundary deformational disp (boundary)
SV = bsxfun(@times,V',S)' ;
%     % Matrix A
% --------
coeffPROY = y'*z ;
[UA,SA,VA] = SVDT(coeffPROY) ;

TOL = 1e-10 ;
if find(SA <1e-10)
    error('This problem is not amenable to reduction. There are "buble" modes ')
end

%IS_ROTATION = 0 ;
%if IS_ROTATION ==1
%    A = UA*VA';
%else
A = UA;
%end

% Coefficients
SVA = SV*A ;
% New matrix


%if INCLUDE_SINGULAR_VALUES ==2

BasisUnew = AsnapDISPdef*SVA ;
% else
%     BasisUnew =  BasisUdef_L*SVA ;
%
% end


if size(BasisUnew,2) ~=size(z,2)
    
    disp(['Check this routine. Number of displacement modes less than number of reactions'])
    disp(['Set a lower tolerance for the SVD of def. displacements'])
    error(' ')
end
%
