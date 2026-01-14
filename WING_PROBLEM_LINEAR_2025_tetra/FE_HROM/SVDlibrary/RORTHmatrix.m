function [Q B GAMMA,alpha_b,beta_b,Rmax,ETIME_row] = RORTHmatrix(A,epsilon,DATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a block matrix A, RORTHmatrix returns a cell array of basis matrices, Q, that approximately
% spans the  column space of A (with accuracy specified by epsilon).
% --------------------------
% INPUTS  (mandatory)
%------------------------------------------------------------------------------
%  A:      p x q cell array containing a conforming partition of A, i.e.:
%
%     A = {A_11, A_12 ... A_1q
%       A_21, A_22 ....A_2q
%       ... .... ..... ...
%       A_p1, A_p2 ... A_pq}
%
%     In turn, each entry A_ij may be a numeric array, or the name of a MAT-file containing
%     the numeric array. The second option is preferred when the whole
%     matrix does not fit into fast memory
%
%   epsilon : Matrix of threshold errors.  size(epsilon) = size(A), and
%     norm(E{i,j},'fro')/norm(A{i,j},'fro') <= epsilon(i,j), if
%     DATA.EPSILON_ABSOLUTE =0 , then
%     norm(E{i,j},'fro') <= epsilon(i,j)
%
%     where E{i,j} = A{i,j} -  Q{i}*Q{i}'*A{i,j}
%
%
%   OPTIONAL INPUTS  (with default values )
%   ----------------
%   1)  DATA.EPSILON_ABSOLUTE   = 0  (see definition epsilon)
%   2)  DATA.Rini = 0
%   3)  DATA.rho_est =0.05
%   4)  DATA.PARALLEL = 0
%
%    DATA.Rini = Initial estimation for the rank of submatrix A{1,1}. The
%    default value is min(DATA.rho_est*ceil(max(size(A{1,1})))).
%
%    DATA.PARALLEL = Enable parallelization of the loop over row blocks.
%    The only difference of the parallel loop and the serial one is that,
%    in the serial one, the maximum rank of the submatrices of A(:,i) is
%    used as initial estimation for the rank of A(1,i+1)
%
%    NOTE: If you try to run with the parallel option, a warning message
%    will appear indicating that the parallelization configuration should
%    be included before the parfor loop (the default options may not be optimal).
%
% -----------------------------------------------------------------------------
%
%  OUTPUTS
%  ------
%  Q  -->   (1 x size(A,1)) cell array such that Q{i}  approximately spans the
%  column space of row block  A(i,:)  (see description input data epsilon)%
%
%   B -->   Cell array of dimensions size(B) = size(A'),  where B{j,i} =
%   A{j,i}'*Q{i}
%
%   GAMMA --> size(A) matrix containing the difference between the error made in approximating
%   each submatrix, and the specified tolerance:
%   GAMMA(i,j) = sqrt(epsilon(i,j)^2*norm(A{i,j}'fro')^2 - norm(E{i,j}),'fro')
%
%   alpha_b and beta_b: Row matrices with the dimensions of B  (size(B{i,j}) = [alpha_b(i) beta_b(j)])
%
%   Rmax: max(rank(A{i,j})). Maximum rank of all submatrices of A
%
%  ETIME_row --> Cell array containing the times required for computing
%  each Q{i}
%
%
%  Written by Joaquín A. Hernández Ortega, Dec. 2016
%  UPC/CIMNE. jhortega@cimne.upc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'EPSILON_ABSOLUTE',0) ;
DATA = DefaultField(DATA,'Rini',0) ;
DATA = DefaultField(DATA,'alpha',[]) ;
DATA = DefaultField(DATA,'beta',[]) ;
DATA = DefaultField(DATA,'normAoriginal',[]) ;
DATA = DefaultField(DATA,'HIDE_OUTPUT',1) ;

  



p = size(A,1) ; Q = cell(1,p) ;   B = cell(size(A')) ; q = size(A,2) ;
GAMMA = zeros(size(A)) ; alpha_b = zeros(1,q) ; beta_b = zeros(1,p) ;
Rmax = 0 ; ETIME_row =cell(1,p);
%dbstop('17')
if DATA.HIDE_OUTPUT == 0
disp('----------------')
disp('LOOP OVER ROW BLOCKS')
disp('----------------')
end

%dbstop('91')
ncolMAX = zeros(1,p) ;
if DATA.PARALLEL == 0
    for i=1:p
        if DATA.HIDE_OUTPUT == 0

        disp('------------------------------------------')
        disp(['ROW block = ',num2str(i), '  of ',num2str(p)])
        disp('------------------------------------------')
        end
        
        if ~isempty(DATA.normAoriginal)
            try
            normB = DATA.normAoriginal(i,:) ;
            catch 
                normB = [] ; 
            end
        else
            normB = [] ; 
        end
        [Q{i},B(:,i),GAMMA(i,:),nrowA,beta_b(i),Rmax_loc,ETIME_row{i}] = RORTH1row(A(i,:),epsilon(i,:),DATA,normB)   ;
        DATA.Rini = Rmax_loc ;
        Rmax = max(Rmax,Rmax_loc) ;
     %   beta_b(i) = size(Q{i},2) ;
    end
    alpha_b = nrowA ;
else
    
 %   error('Set below the parallelization parameters (or comment this line  if you want to use the default ones)')
    if matlabpool('size') ==0
        matlabpool open
    end
  %  matlabpool open
    nrowA = cell(p,1); Rmax_loc = zeros(p,1) ;  
    parfor i=1:p
        
          if isempty(DATA.normAoriginal)
            normB = DATA.normAoriginal(i,:) ;
        else
            normB = [] ; 
        end
        
        [Q{i},B(:,i),GAMMA(i,:),nrowA{i},beta_b(i),Rmax_loc(i),ETIME_row{i}] = RORTH1row(A(i,:),epsilon(i,:),DATA,normB)   ;
     %   beta_b(i) = size(Q{i},2) ;
        
    end
    
 %   matlabpool close
    alpha_b = nrowA{1} ;
    Rmax = max(Rmax_loc) ;
end

end