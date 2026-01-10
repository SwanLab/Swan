function      Q  =RodriguesFormula(a)  
% Rodrigues formula, rotation
% - JAHO, 25-FEB-2025, UPC,TERRASSA 
if  nargin == 0
    a = [0,0,0.1]' ; 
end

% % 
%  \begin{equation}
%  \label{eq:::3}
% %     Q =    \ident + \seno{\normd{a}} \spin{a/\normd{a}} + (1- \coseno{\normd{a}})  \spin{a/\normd{a}}^2
%  \end{equation}

% Norm of a 
ndim  =3;  
nelem  = length(a)/ndim ; 

% Computing the normalized angle 
a_matrix = reshape(a,ndim,[]) ; 
% Norm of a 
normd_a = sqrt(sum(a_matrix.^2,1));

a_norm = bsxfun(@times,a_matrix',1./normd_a')' ; 

% SPIN
S_anorm = SpinVector(a_norm) ;  

% IDENTITY CONTRIBUTION 
% Qident =    \ident 
Qident = repmat(eye(3),nelem,1) ; 
% SINE CONTRIBUTION
% Qseno =     \seno{\normd{a}} \spin{a/\normd{a}}
term_seno_normd_a = sin(normd_a) ; 
term_seno_normd_a = repmat(term_seno_normd_a,ndim,1) ; 
term_seno_normd_a = term_seno_normd_a(:) ; 
Qseno = bsxfun(@times,S_anorm,term_seno_normd_a) ; 
% COSINE CONTRIBUTION
% Qcoseno =  (1- \coseno{\normd{a}})  \spin{a/\normd{a}}^2
term_coseno_normd_a = 1-cos(normd_a) ; 
term_coseno_normd_a = repmat(term_coseno_normd_a,ndim,1) ; 
term_coseno_normd_a = term_coseno_normd_a(:) ; 
S_anorm_2 = MultiplyMatrixBlocks(S_anorm,S_anorm) ; 
Qcoseno = bsxfun(@times,S_anorm_2,term_coseno_normd_a) ;
Q = Qident + Qseno + Qcoseno ; 

