
function [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(BasisF,SingVal_F,W,DATA)
% Version before 4th Dec. 2019

if nargin == 0
    %DATA = [];
    load('tmp1.mat')
    %     DATA.IncludeSingularValuesF =0 ;
    %     DATA.TOL = 1e-3 ;
    % G = J ;
end

if isempty(SingVal_F) || length(SingVal_F) ~= size(BasisF,1)
    SingVal_F  =ones(size(BasisF,2),1) ;
end

DATA = DefaultField(DATA,'IncludeSingularValuesF',0) ; %  ; .IncludeSingularValuesF = 1
if DATA.IncludeSingularValuesF == 1
    warning('This option has proved unreliable...disable it')
    %   G = bsxfun(@times,BasisF', (SingVal_F));  %  Version before  5th Dec-2019...
    G = bsxfun(@times,BasisF',sqrt(SingVal_F));  %
    b = G*sqrt(W) ;  % b Vector (exact integral)
    bEXACT = b ;
else
    G = BasisF' ;
    b = G*sqrt(W) ;  % b Vector (exact integral)
    if ~isempty(SingVal_F)
        bEXACT = b.*SingVal_F ;
    else
        bEXACT =b ;
    end
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



DATA = DefaultField(DATA,'IND_POINTS_CANDIDATES',[]) ;

if ~isempty(DATA.IND_POINTS_CANDIDATES)
    y = intersect(y,DATA.IND_POINTS_CANDIDATES) ;
end
yORIG = y ; 

DATA  = DefaultField(DATA,'USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR',1) ;  % 27th-april-2020
USEsingvERR = DATA.USE_SINGULAR_VALUES_INTF_FOR_MEASURING_ECM_ERROR ; 


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
y = y(:) ;
NITERACIONES = 10*m ; 
 
while  nerrorACTUAL >TOL && mPOS <m   && ~isempty(y) && k<NITERACIONES
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
    if  ~isempty(n)
        % STEP 5
        y = [y; z(n)];  z(n)=[] ;
        H = MultiUpdateInverseHermitian(H,n) ;
        % Recomputing alpha
        alpha = H*(G(:,z)'*b );
    end
    % STEP 6
    r = b-G(:,z)*alpha ;
    nerror = norm(r)/norm(b) ; % Relative error (using r and b)
    if DATA.IncludeSingularValuesF == 0 && USEsingvERR == 1
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
    
    %     if length(z) == m
    %         dbstop('88')
    %         disp('')
    %     end
    
end


if  k>= NITERACIONES
    PROPORpoints = length(yORIG)/length(W)*100; 
    error(['NO CONVERGE. ENLARGE THE SET OF CANDIDATE POINTS  (NOW IT CONTAINS ',num2str(PROPORpoints),' % of the total number of Gauss points)'])
end 


w = alpha.*sqrt(W(z)) ;

disp(['Total number of iterations =',num2str(k)])

figure(500)
hold on
xlabel('Number of points')
ylabel('Error (%)')
plot(NPOINTS,ERROR_GLO*100,'k')

DATAOUT.kiteraciones = k ;

% figure(501)
% hold on
% xlabel('Number of iterations')
% ylabel('Number of points')
% plot(NPOINTS,'k')


