function [z,w]= EmpiricalCubatureMethod_NNLS(BasisF,SingVal_F,W,DATA)
% Version after  4th Dec. 2019

if nargin == 0
    %DATA = [];
    load('tmp.mat')
    DATA.IncludeSingularValuesF =1 ; 
    DATA.TOL = 1e-6 ; 
    % G = J ;
end

DATA = DefaultField(DATA,'IncludeSingularValuesF',1) ; %  ; .IncludeSingularValuesF = 1
if DATA.IncludeSingularValuesF == 1
 %   G = bsxfun(@times,BasisF', (SingVal_F));  %  Version before  5th Dec-2019...  
    G = bsxfun(@times,BasisF',sqrt(SingVal_F));  % 
    b = G*sqrt(W) ;  % b Vector (exact integral)
    bEXACT = b ;
else
    G = BasisF' ;
    b = G*sqrt(W) ;  % b Vector (exact integral)
    bEXACT = b.*SingVal_F ; 
end
    nbEXACT = norm(bEXACT) ; 


Gnorm =sqrt(sum(G.*G,1)) ; % Norm of Columns of G
M = size(G,2) ;  % Number of FE points
DATA = DefaultField(DATA,'TOL',0) ; % Default tolerance for convergence
TOL = DATA.TOL ;
% INITIALIZATIONS
% ------------------------
z = [] ; % Set of integration points
% Set of candidate points (those whose associated column has low norm are removed)

 %  PointsWithZero =  find(sum(G(1:end-1,:),1)==0) ; 

y=1:M ;
%y(PointsWithZero) = []  ; 
DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-6) ;
GnormNOONE =sqrt(sum(G(1:end-1,:).*G(1:end-1,:),1)) ; % Norm of Columns of G

if DATA.TOLFilterCandidatePoints >0
    TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
    rmvpin = find(GnormNOONE(y)<TOL_REMOVE) ;
    y(rmvpin) = [] ;
end

DATA = DefaultField(DATA,'RemoveColumnsWithNegativeProjection',0); % 4-Dec-2019


 






alpha = [] ; % Vector of weights
mPOS = 0 ; % Number of nonzero weights
r = b ; % Residual vector
k = 1;  % Number of iterations
errorGLO = [] ; %(for storing error)
% Default number of points
DATA = DefaultField(DATA,'npoints',length(b)) ;
m = min(DATA.npoints,length(b)) ;
% END INITIALIZATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normB = norm(b) ;
nerror = norm(r)/normB  ;
H = [] ; % Inverse of (Gz'*Gz)
ERROR_GLO = [] ;
NPOINTS =[] ; 
nerrorACTUAL = nerror; 
WEIGHTS_old = zeros(size(G,2),1); 
WEIGHTS_new = zeros(size(G,2),1); 

 
while  nerrorACTUAL >TOL && mPOS <m   && ~isempty(y)
    % TIMELOC_k = tic ;
    % STEP 1. Compute new point
    ObjFun = G(:,y)'*r ;
    ObjFun = ObjFun./Gnorm(y)';
    [maxLOC, indSORT] = max(ObjFun)  ;
    i = y(indSORT(1)) ;
    % STEP 2.  Update alpha and H  (unrestricted least-squares)
    if k==1
        alpha =  G(:,i)\b ;
        H = 1/(G(:,i)'*G(:,i)) ;
    else
        [H alpha] = UpdateWeightsInverse(G(:,z),H,G(:,i),alpha,r) ;
    end
    % STEP 3. Move i from set y to set z
    z = [z;i] ;     y(indSORT(1)) = [] ;
    % STEP 4. Find possible negative weights
    n = find(alpha<=0) ;
   
     WEIGHTS_new(z) = alpha ; 
     ineg = 1; 
    while  ~isempty(n)
        % STEP 5. Replaced by a while     
       
      %  eta  = WEIGHTS_old(z(n))./(WEIGHTS_old(z(n))-WEIGHTS_new(z(n))) ;  
      %  eta = min(eta) ; 
      %  WEIGHTS_old =  WEIGHTS_old + eta*(WEIGHTS_new-WEIGHTS_old) ; 
      if length(n) >1
          disp('')
      end
        [aaaa bbbbb]=  sort(alpha(n)) ;  % Sort in ascending order
        n = n(bbbbb);
        
        y = [y z(n(ineg))'];  z(n(ineg))=[] ;
        H = MultiUpdateInverseHermitian(H,n(ineg)) ;
        % Recomputing alpha
        alpha = H*(G(:,z)'*b );
             n = find(alpha<=0) ;
        
    end
    % STEP 6
    r = b-G(:,z)*alpha ;
    nerror = norm(r)/norm(b) ; % Relative error (using r and b)
    if DATA.IncludeSingularValuesF == 0
        nerrorACTUAL = SingVal_F.*r ;
        nerrorACTUAL = norm(nerrorACTUAL/nbEXACT); 
    else
        nerrorACTUAL = nerror ;
    end
    % STEP 7
    disp(['k = ',num2str(k),', m=',num2str(length(z)),' ,','; error n(res)/n(b) (%) = ',...
        num2str(nerror*100),';  Actual error % =',num2str(nerrorACTUAL*100)]) ; 
    ERROR_GLO(k) = nerrorACTUAL ;
    NPOINTS(k) =  length(z) ; 

    mPOS = length(z) ; 
    k = k + 1 ;
    WEIGHTS_old = WEIGHTS_new ; 
    
%     if length(z) == m
%         dbstop('88')
%         disp('')
%     end
   
end





w = alpha.*sqrt(W(z)) ;

disp(['Total number of iterations =',num2str(k)])

figure(500)
hold on
xlabel('Number of points')
ylabel('Error (%)')
plot(NPOINTS,ERROR_GLO*100,'k')

% figure(501)
% hold on
% xlabel('Number of iterations')
% ylabel('Number of points')
% plot(NPOINTS,'k')


