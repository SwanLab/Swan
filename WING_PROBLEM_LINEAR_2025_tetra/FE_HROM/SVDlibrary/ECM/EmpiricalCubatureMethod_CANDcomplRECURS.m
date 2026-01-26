function [z,w,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcomplRECURS(BasisF,SingVal_F,W,DATA)
% Version created 6-OCT-2021
% Similar to EmpiricalCubatureMethod_CANDcompl.m (complementary methods)
% More details can be found in
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/OptimalECM/PolyIMPLEMENT.m
% ----------------------------------------------------
error('This option hasn not proved to offer any advantage ....')
if nargin == 0
    
    load('tmp.mat')
    
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
DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-10) ;
GnormNOONE =sqrt(sum(G(1:end-1,:).*G(1:end-1,:),1)) ; % Norm of Columns of G

if DATA.TOLFilterCandidatePoints >0
    TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
    rmvpin = find(GnormNOONE(y)<TOL_REMOVE) ;
    y(rmvpin) = [] ;
end

DATA = DefaultField(DATA,'RemoveColumnsWithNegativeProjection',0); % 4-Dec-2019



DATA = DefaultField(DATA,'IND_POINTS_CANDIDATES',[]) ;

if ~isempty(DATA.IND_POINTS_CANDIDATES)
    idomain = 1;
    % We first try with just the fist domain
    y = intersect(y,DATA.IND_POINTS_CANDIDATES{idomain}) ;
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
iOLD = [] ;
%yCOMPL = (1:length(W))';
%yCOMPL(y) = [] ; % Complementary set of the candidate set (in case this set is not sufficiently rich)


NITERATIONS_NO_MORE_POINTS = 10 ; %    NUMBER OF ITERATIONS IN WHICH THE CODE IS ALLOWED TO "ITERATE" IN ORDER TO FIND
%  THE MINIMUM

ITER_MINIMUM = 0;
MAX_NUMBER_POINTS = 0 ;

yCOMPL = [] ;


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
        if ITER_MINIMUM > NITERATIONS_NO_MORE_POINTS
            disp('--------------------------------------------------------------')
            disp(['The algorithm cannot proceed with the current set of points (up to domain =',num2str(idomain),')'])
            
            
            if  idomain < length(DATA.IND_POINTS_CANDIDATES)
                disp(['We enlarge the set of candidates with the complementary set up to domain =',num2str(idomain+1)])
                yCOMPL  = [yCOMPL; DATA.IND_POINTS_CANDIDATES{idomain+1}'] ;
                y = [y;yCOMPL] ;
                idomain = idomain + 1;
            end
        end
        % ITER_MINIMUM = ITER_MINIMUM + 1;
    else
        % ITER_MINIMUM = 0 ;
    end
    
    if length(z) > MAX_NUMBER_POINTS
        ITER_MINIMUM = 0 ;
    else
        ITER_MINIMUM = ITER_MINIMUM + 1;
    end
    
    iOLD = i ;
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
    MAX_NUMBER_POINTS = max(mPOS,MAX_NUMBER_POINTS) ;
    k = k + 1 ;
    
    %     if length(z) == m
    %         dbstop('88')
    %         disp('')
    %     end
    
end




if  k>= NITERACIONES
    PROPORpoints = length(yORIG)/length(W)*100;
    error(['NO CONVERGENCE. ENLARGE THE SET OF CANDIDATE POINTS  (NOW IT CONTAINS ',num2str(PROPORpoints),' % of the total number of Gauss points)'])
end


w = alpha.*sqrt(W(z)) ;

disp(['Total number of iterations =',num2str(k)])

PLOT_error = 0 ;

if PLOT_error == 1
    
    figure(500)
    hold on
    xlabel('Number of points')
    ylabel('Error (%)')
    plot(NPOINTS,ERROR_GLO*100,'k')
end

DATAOUT.kiteraciones = k ;

% figure(501)
% hold on
% xlabel('Number of iterations')
% ylabel('Number of points')
% plot(NPOINTS,'k')


