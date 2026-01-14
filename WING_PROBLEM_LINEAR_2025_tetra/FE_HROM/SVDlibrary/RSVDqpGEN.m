function [U,S,V,ETIME,RankR,e_svd] = RSVDqpGEN(A,epsilon,DATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a matrix A, RSVDqpGEN returns   a factorization  [U,S,V] of the form
%   A =  U*S*V^T + E
% --------------------------
% INPUTS
%------------------------------------------------------------------------------
%  A:    Matrix to be factorized. It can be given as a single numeric array,
%       or as a  p x q cell array containing a conforming partition of A, i.e.:
%
%     A = {A_11, A_12 ... A_1q
%       A_21, A_22 ....A_2q
%       ... .... ..... ...
%       A_p1, A_p2 ... A_pq}
%
%     In turn, each entry A_ij may be a numeric array,  the name of a MAT-file containing
%     the numeric array or a cell array of  mat-files.  
%
%   OBSERVATION: In the latter case, when any of the entries of A is a cell
%   array of mat files, the two-level function  RSVDqpTWOL is invoked.
%   Otherwise, function RSVDqp is employed.  
%
% -------------------------------------------------------------------------------
%  epsilon :
%     If A is a single numeric array, epsilon is a positive number (below 1) that
%     indicates the error threshold in approximating the whole matrix
%     (in the frobenius norm), i.e.
%     norm(E,'fro')/norm(A,'fro') <= epsilon, where E = A - U*S*V^T
%
%     If A is a cell array, size(epsilon) = size(A), and
%     norm(E{i,j},'fro')/norm(A{i,j},'fro') <= epsilon(i,j)
%
%     NOTES: 
%     a) If all the entries of epsilon are zero, then the global error has to
%     be specified by variable  DATA.EPSILON_GLO   (by default = 0)
%     
%     b) If DATA.EPSILON_GLO_isRELATIVE =0, the global threshold is given
%     in absolute value
%
% -----------------------------------------------------------------------------
%
%  OUTPUTS
%  ------
%  U --> (M x r)  Matrix of left singular vectors (it approximately spans
%  the column space of A); r denotes the rank of the approximation.
%  V --> (N x r)  Matrix of right singular vectors (it approximately spans
%  the row space of A)
%  S ---> Vector of singular values  (r x 1)
%  ETIME --> Structure containing the times required for each operation
%
%  NOTES: If epsilon(i,j) = 0 , then the obtained
%  factorization is the classical truncated SVD
%  ----------------
%  OPTIONAL INPUTS, (with default values)
%  -----------------
% 1)  DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP = 0
% 2)  DATA.Rini = 0 
% 3)  DATA.rho_est = 0.05 
% 4)  DATA.TRANSPOSE = 0    (only employed in the  one-level method)
% 5)  DATA.COMPUTE_V_LAST_SVD = 1
% 6)  DATA.TWOLEVELMETHOD.LimitGbytes =3; 
% 
%
%  When each A{i,j} is a .MAT
%  filename,    the contents of A{i,j} are discarded from memory
%  after its corresponding incremental basis vectors are computed.
%  The product of the transponse of basis matrix Qi  times the i-th row
%  block is stored in memory if variable
%  DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP = 1 ;
% 
%    DATA.Rini = Initial estimation for the rank of submatrix A{1,1}. The
%    default value is min(DATA.rho_est*ceil(max(size(A{1,1})))).   
%
%  DATA.TRANSPOSE -> The factorization is calculated on the transpose of
%  the matrix(only employed in the  one-level method)
%
%  DATA.COMPUTE_V_LAST_SVD -->  If 0, the matrix of right singular vectors
%  is not calculated
%
%  DATA.TWOLEVELMETHOD.LimitGbytes =3. When the two-level method is
%  invoked, this variable indicates the maximum size (in Gbytes) of a submatrix A{i,j}
%  for being permanently loaded in the RAM  
%  
%  
%
%  Written by Joaquín A. Hernández Ortega, Dec. 2016, January 2017
%  UPC/CIMNE. jhortega@cimne.upc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('DefaultField') == 0
    addpath('AUXFUN')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT INPUTS
% ----------------------------------------------------
if prod(size(epsilon)) == 1
DATA = DefaultField(DATA,'EPSILON_GLO',epsilon) ; % See preabmble
else
    DATA = DefaultField(DATA,'EPSILON_GLO',0) ;
end
DATA = DefaultField(DATA,'USE_SLOW_MEMORY',[]) ;  %

DATA.USE_SLOW_MEMORY = DefaultField(DATA.USE_SLOW_MEMORY,'STORE_ALSO_SECOND_STEP',0) ;  %
DATA = DefaultField(DATA,'EPSILON_GLO_isRELATIVE',1) ; % For multi-level version
DATA = DefaultField(DATA,'COMPUTE_V_LAST_SVD',1) ; % It computes the matrix of right-singular vectors
DATA = DefaultField(DATA,'Rini',0) ; % Initial estimation for the rank
DATA = DefaultField(DATA,'TRANSPOSE',0) ; % 
DATA = DefaultField(DATA,'PARALLEL',0) ; % 
DATA = DefaultField(DATA,'rho_est',0.05) ; % 

 DATA = DefaultField(DATA,'COMPLETE_SVD',0) ; 


DATA = DefaultField(DATA,'TWOLEVELMETHOD',[]) ; 
DATA.TWOLEVELMETHOD = DefaultField(DATA.TWOLEVELMETHOD,'LimitGbytes',3) ; 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 

% Checkjwhether the two-level method comes into play
if iscell(A)
     CHECK = cellfun(@iscell,A) ;
     if  any(CHECK)
         DATA.TWOLEVELMETHOD.ACTIVE = 1; 
     else
         DATA.TWOLEVELMETHOD.ACTIVE = 0; 
     end    
end
RankR = 0 ; 
if DATA.TWOLEVELMETHOD.ACTIVE == 0
    % One-level method
    if DATA.COMPLETE_SVD == 0 | all(epsilon==0)     
    [U,S,V,ETIME,e_svd,RankMatrix] = RSVDqp(A,epsilon,DATA) ;
    else
         [U,S,V,ETIME,RankR,e_svd] = RSVDqpCOMP(A,epsilon,DATA) ;
    end
    
    
else
  %  Two-level method
    [U,S,V,ETIME] = RSVDqpTWOL(A,epsilon,DATA) ;
end



 
