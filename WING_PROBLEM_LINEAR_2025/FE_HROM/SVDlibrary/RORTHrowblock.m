function [Q,L,gamma,NormB,ETIME] = RORTHrowblock(B,epsilon,DATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a row block matrix B, RORTHrowblock returns a matrix Q  that approximately
% spans the  column space of A (with accuracy specified by epsilon).
% --------------------------
% INPUTS  (mandatory)
%------------------------------------------------------------------------------
%  B:      1 x q cell array containing a conforming partition of B, i.e.:
%
%     B = {B_11, B_12 ... B_1q}
%
%     In turn, each entry B_ij may be a numeric array, or the name of a MAT-file containing
%     the numeric array. The second option is preferred when the whole
%     matrix does not fit into fast memory
%
%   epsilon : Row matrix of threshold errors.  size(epsilon) = size(B), and
%     norm(E{i},'fro')/norm(B{i},'fro') <= epsilon(i). If
%     DATA.EPSILON_ABSOLUTE =0 , then
%     norm(E{i},'fro') <= epsilon(i)
%
%     where E{i} = B{i} -  Q*Q'*B{i}
%
%
%   OPTIONAL INPUTS  (with default values )
%   ----------------
%   1)  DATA.EPSILON_ABSOLUTE   = 0  (see definition epsilon)
%   2)  DATA.Rini = 0
%   3)  DATA.rho_est =0.05
%
%    DATA.Rini = Initial estimation for the rank of submatrix B{1}. The
%    default value is min(DATA.rho_est*ceil(max(size(B{1})))).
%
%
% -----------------------------------------------------------------------------
%
%  OUTPUTS
%  ------
%  Q  -->   Orthogonal matrix that   approximately spans the
%  column space of     B  (see description input data epsilon)%
%
%   GAMMA --> size(B) matrix containing the difference between the error made in approximating
%   each submatrix, and the specified tolerance:
%   GAMMA(i) = sqrt(epsilon(i)^2*norm(B{i}'fro')^2 - norm(E{i}),'fro')
%
%   nrowL: Number of rows of each submtrix  of L
%
%
%  ETIME --> Structure array containing the times required for computing Q
%
%
%  Adaptation of RORTH1row.m (in Oct-2022)
%  See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Fun3D/IssueSVD.mlx
%  UPC/CIMNE. jhortega@cimne.upc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('59')
if nargin == 0
    load('tmp.mat')
end

normB = [] ;

ETIME.loop = tic ;
DATA = DefaultField(DATA,'EPSILON_ABSOLUTE',0) ;
DATA = DefaultField(DATA,'Rini',0) ;  % Estimation for the initial rank
DATA = DefaultField(DATA,'rho_est',0.05) ;  % Parameter  0<rho<1 for estimating the initial rank
%DATA = DefaultField(DATA,'TryToPrecomputeQtB',1) ;
%DATA = DefaultField(DATA,'NoCalculateQTAWhenEPS1',0) ;
DATA = DefaultField(DATA,'HIDE_OUTPUT',1) ;



%---------------------------------------------------------------------------------------
R = zeros(size(B)) ;
q = length(B) ;  % Number of submatrices
Rest = DATA.Rini ; % Estimation for the rank of the first matrix to be treated
% Submatrices with lower tolerances are addressed first
[~, INDICES] = sort(epsilon) ;
epsilonSVD = zeros(size(INDICES)) ; % Array with approximation errors
gamma =  zeros(size(INDICES)) ;
%bnorm = gamma ;
Q = [] ;  % Initialization orthogonal matrix
L = cell(1,q) ; % Cell array containing Q'*B{i}
% ncolMAX  = 0;
% nrowL  = zeros(size(B)) ;
NormB = zeros(size(B)) ; 
 


if DATA.HIDE_OUTPUT == 0
    disp('----------------')
    disp('LOOP OVER COLUMN SUBMATRICES')
    disp('----------------')
end
for j=1:q
    i = INDICES(j) ; % We begin with the submatrix with lower tolerance
    if DATA.HIDE_OUTPUT == 0
        disp('------------------------------------------')
        disp(['i = ',num2str(j), '  of ',num2str(q)])
    end
    
    if isempty(B{i}) | ((epsilon(i) == 1) && DATA.SetGammaToZero ==1 )
        if DATA.HIDE_OUTPUT == 0
            disp('Ignoring matrix')
        end
    else
        
        if ischar(B{i})
            % We have to load from memory the corresponding numeric array
            if DATA.HIDE_OUTPUT == 0
                disp(['Loading ...'])
            end
            SSS = load(B{i}) ; fff = fieldnames(SSS) ;
            %        if DATA.TRANSPOSE == 0
          %  Bi = SSS.(fff{1}) ; 
            %       else
            %          Bi = SSS.(fff{1})' ; SSS = [] ;
            %      end
            
            if ~isempty(DATA.WeightsPremultipy_matrix)
                 Bi = bsxfun(@times,SSS.(fff{1}),DATA.WeightsPremultipy_matrix) ;
            else
                Bi = SSS.(fff{1}) ;
            end
            SSS = [] ;
            if DATA.HIDE_OUTPUT == 0
                disp(['...Done'])
            end
        else
            %  if DATA.TRANSPOSE == 0
          
             if ~isempty(DATA.WeightsPremultipy_matrix)
                 Bi = bsxfun(@times, B{i} ,DATA.WeightsPremultipy_matrix) ;
             else
                   Bi = B{i} ;
            end
            %     else
            %        Bi = B{i}' ;
            %   end
        end
        %%% Estimating the rank of the matrix
        Ni = size(Bi,2) ;     Mi = size(Bi,1) ;
        RiUP = min(Ni,Mi); % Upper bound for the rank
        if Rest ==0; Rest = ceil(DATA.rho_est*RiUP);   end
        Rest = min(Rest,RiUP);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j==1 || isempty(Q)
            dB = Bi ;
        else
            if DATA.HIDE_OUTPUT == 0
                disp('Computing residual ...dB = Bi - Q*(Q^T*Bi)')
            end
            dB = Bi - Q*(Q'*Bi) ;
            if DATA.HIDE_OUTPUT == 0
                disp('...done')
            end
        end
        
        %   if isempty(normB)
        b = norm(Bi,'fro') ;
        NormB(i) = b; 
        %  else
        %     b= normB(i) ;
        % end
        mu = (max(size(dB))*eps(b))  ;% Machine precision (calculated with the original matrix, not the residual)
        if epsilon(i) > 0
            if DATA.EPSILON_ABSOLUTE ==0
                e0 = epsilon(i)*b ;
            else
                e0 = epsilon(i) ;
            end
        else
            % If e0, then we set the absolute tolerance to the machine
            % precision parameter
            e0 = mu ;
        end
        % If e0 is not zero, but happens to be smaller than the machine
        % precision parameter, we then set e0 = max(e0,mu)
        e0 = max(e0,mu) ;
        
        if ~isempty(dB)
            if j>1
                db = norm(dB,'fro'); % Norm of the residual
            else
                db = b ;
            end
            
            if db > e0 || j==1
                % If the norm of the residual is larger than the prescribed
                % tolerance, then we iteratively compute the new basis vectors
                % of Q by invoking function RSVDt
                DATA.Rest  =Rest ;
                DATA.COMPUTE_V_SVD = 0 ;
                DATA.RELATIVE_SVD = 0;   % Modification 4-March-2020 !
                [dQ,Hi,~,eSVD,RankMatrix] = RSVDT(dB,e0,mu,Rest,DATA);
                if isempty(dQ)
                    eSVD = db ;
                end
            else
                dQ = [] ; Hi = [] ;
                eSVD = db ;
                RankMatrix = Rest ;
            end
            
            if epsilon(i) > 0
                epsilonSVD(i) = eSVD/b ;
                gamma(i) =  sqrt((epsilon(i)*b)^2 - eSVD^2) ;
            else
                epsilonSVD(i) = eSVD ;
            end
        else
            dQ = [] ;        Hi = [] ;                eSVD = [] ;        RankMatrix = Rest ;
        end
        
        if j>1 && ~isempty(dQ) && ~isempty(Q)
            % Re-orthogonalization
            dQ = orth(dQ - Q*(Q'*dQ)) ;
        end
        % Augmenting the basis matrix
        Q = [Q dQ];
        
        % Estimation for the next submatrix
        if   ~isempty(dQ)
            Rest = RankMatrix ;
        else
            Rest = 1 ;
        end
        
        if DATA.HIDE_OUTPUT == 0
            disp(['K = ',num2str(length(Hi)),' of ',num2str(Ni),' columns (rank = ',num2str(RankMatrix),')'])
            disp(['ERROR = ', num2str(epsilonSVD(i)),' (THRESHOLD:',num2str(epsilon(i)),' )'])
        end
        R(j) = RankMatrix;
        
        %%%
        if  epsilon(i)==0
            % Pre-computing L
            %     if isstr(Bi)
            %                 if  DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP ==0
            %                     L{i} =Q'*Bi ;
            %                     nrowL(i) = size(L{i},1) ;
            %                     ncolMAX = max(ncolMAX,size(L{i},2)) ;
            %                 else
            %                     dirREF = dir(Bi) ; dirREF = dirREF.name ;
            %                     L{i}  =  [DATA.PATH_STORE_AUX,dirREF(1:end-4),num2str(i),'.mat'] ;
            %                     Li = Bi'*Q ;
            %                     save(L{i},'Li') ;
            %                     nrowL(i) = size(Li,1) ;
            %                     ncolMAX = max(ncolMAX,size(Li,2)) ;
            %                 end
            %             else
            %
            %  if  ~isempty(Bi) & ~isempty(Q)
            % if DATA.TRANSPOSE == 0
            L{i} =Q'*Bi;
            %    else
            %       L{i} =Bi*Q ;
            %  end
            % end
            
            
            %    nrowL(i) = size(L{i},1) ;
            %    ncolMAX = max(ncolMAX,size(L{i},2)) ;
            %    end
            
        end
        
    end
    
    %%%
    
end


ETIME.loop = toc(ETIME.loop) ;


ETIME.prodBQ = tic  ;

if DATA.HIDE_OUTPUT == 0
    disp('---------------------------------')
    disp('Matrix L = B^T*Q ...')
    disp('---------------------------------')
end

for i = 1:length(B)
    if DATA.HIDE_OUTPUT == 0
        disp(['i=',num2str(i)])
    end
    if  ~isempty(Q) &&    epsilon(i)~=0
        if isstr(B{i})
            loct1 =  tic ;
            if  DATA.HIDE_OUTPUT == 0
                disp(['Loading ...'])
            end
            SSS = load(B{i}) ;
            fff = fieldnames(SSS) ;
            %    if DATA.TRANSPOSE ==0
            Bi = SSS.(fff{1}) ; SSS = [] ;
            %   else
            %      Bi = SSS.(fff{1})' ; SSS = [] ;
            %  end
            if DATA.HIDE_OUTPUT == 0
                disp(['...Done'])
            end
            %         if  DATA.USE_SLOW_MEMORY.STORE_ALSO_SECOND_STEP ==0
            L{i} =Q'*Bi; % [Q'*B{i}];
            %     ncolMAX = max(ncolMAX,size(L{i},2)) ;
            %    nrowL(i) = size(L{i},1) ;
            %             else
            %                 %  dbstop('173')
            %                 dirREF = dir(B{i}) ; dirREF = dirREF.name ;
            %                 L{i}  =  [DATA.PATH_STORE_AUX,dirREF(1:end-4),num2str(i),'.mat'] ;
            %                 Li = B{i}'*Q ;
            %                 ncolMAX = max(ncolMAX,size(Li,2)) ;
            %                 nrowL(i) = size(Li,1) ;
            %                 save(L{i},'Li') ;
            %             end
        else
            %   if DATA.TRANSPOSE == 0
            L{i} =Q'*B{i} ;
            %  else
            %     L{i} =B{i}*Q ;
            % end
            %ncolMAX = max(ncolMAX,size(L{i},2)) ;
            
            %nrowL(i) = size(L{i},1) ;
        end
        
    end
    
    
end



for i = 1:length(L)
    if size(L{i},1) ~= size(Q,2)
        addZEROS =  size(Q,2)-size(L{i},1) ;
        L{i } = [L{i}; zeros(addZEROS,size(L{i},2))  ] ;
    end
%     if size(L{i},1) ~= DATA.beta(i)
%         addZEROS = DATA.beta(i)-size(L{i},1) ;
%         L{i } = [L{i};  zeros(addZEROS,size(L{i},2))  ] ;
%     end
end

%end

% gamma = zeros(size(gamma)) ;

%
%
% RECALCULATE_GAMMA =0 ;
%
%
% if RECALCULATE_GAMMA ==  1
%
% warning('Prueba')
% newgamma = zeros(size(gamma)) ;
% for i = 1:length(L)
%     if epsilon(i) >0
%    proy = B{i} - Q*Q'*B{i} ;
%    actualerror = norm(proy,'fro') ;
%    b = norm(B{i});
%
%       actualerrorREL = actualerror/b;
%
%
%    newgamma(i) = sqrt((epsilon(i)*b)^2-actualerror^2) ;
%
%     end
%
%
%
% end
%
%
% warning('End Prueba')
%
%
% end








if DATA.HIDE_OUTPUT == 0
    disp('...Done')
end



ETIME.prodBQ = toc(ETIME.prodBQ )  ;

