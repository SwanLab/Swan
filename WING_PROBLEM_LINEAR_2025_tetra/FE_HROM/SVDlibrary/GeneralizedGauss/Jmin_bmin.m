function [J,b,Jnorm,DATAOUT,DATA] = Jmin_bmin(Xf,W,DATA)
% Construction of matrix J and vector b
% --------------------------------------
%dbstop('5')
if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'includeSingularValuesinJ',0); % New option
DATA = DefaultField(DATA,'RANDOMIZED_SVD',1); % New option
DATA = DefaultField(DATA,'rankESTIMATE',0); % New option
DATA = DefaultField(DATA,'PLOTsingVAL_Xf_All',0); % New option
DATA = DefaultField(DATA,'IncludeWeights_Jmin',1); % New option
DATA = DefaultField(DATA,'bMIN_withONE',2); % New option




DATA = DefaultField(DATA,'epsilon',[]); % New option

nsnap = size(Xf,2);  % % Number of snashots
M = size(Xf,1) ;     % Number of integration points
%dbstop('15')


if ~iscell(Xf)
    INTexac = Xf'*W ; % High-fidelity integral of all snapshots
else
    INTexac = cell(size(Xf)) ;
    for imat = 1:length(Xf)
        if  ischar(Xf{1})
            load(Xf{imat},'Aloc') ;
            INTexac{imat} =Aloc'*W ;
        else
            INTexac{imat} =Xf{imat}'*W ;
        end
    end
    
    INTexac = cell2mat(INTexac') ;
end




V = sum(W) ;    % Volume of the domain

% Matrix of zero-integral snapshots (step 3 in Box 5.1)
%if  DATA.SVD_BEFORE_zeroaverage ==0
if DATA.IncludeWeights_Jmin ==1
    if DATA.ImposeVolumeConstraint == 1
        Xf = bsxfun(@minus,Xf',INTexac/V)' ;
        Xf = bsxfun(@times,Xf,sqrt(W)) ;
    else
        if ~iscell(Xf)
            Xf = bsxfun(@times,Xf,sqrt(W)) ;
        else
            
            for ipart = 1:length(Xf)
                Xf{ipart} = bsxfun(@times,Xf{ipart},sqrt(W)) ;
            end
            
        end
    end
    
end
%elseif   DATA.SVD_BEFORE_zeroaverage ==1
%   Xf = bsxfun(@times,Xf,sqrt(W)) ;
%  if DATA.MULTIPLY_BY_INTEGRAL == 1
%     INTexacVOL = abs(INTexac) + 1 ;
%     Xf = bsxfun(@times,Xf',INTexacVOL)' ;
% end
%else
%   error('Not implemented option')
%end

% Orthogonal basis matrix Lambda for the column space of Xf (step 4 in Box 5.1)
disp('SVD of Xf ')
%dbstop('20')
%if  DATA.NOsvdXf == 0

% NEWCODE = 1 ;

%dbstop()
%if NEWCODE == 1
%dbstop('73')
if DATA.RANDOMIZED_SVD == 0
    [Lambda,S,V,ERRORsvd,TIMETOTAL] = SVD_BLOCKWISE(Xf,DATA);
else
    DATA.EPSILON_GLO = DATA.TOLsvdXf ;
    if  ~iscell(Xf)
        %  [Lambda,S,V,ERRORsvd,TIMETOTAL,DATA] = SVD_BLOCKWISE_random(Xf,DATA);
        
        [Lambda,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    RSVDqp(Xf, DATA.EPSILON_GLO ) ;
    else
        %     dbstop('80')
        
        
        if isempty(DATA.epsilon)
            epsilon = 0*ones(size(Xf)) ;
        else
            epsilon = DATA.epsilon ;
        end
        [Lambda,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    RSVDqp(Xf,epsilon,DATA) ;
    end
    
    
    
end

DATA = DefaultField(DATA,'TruncationCriterionSVDmodesJumpOrderMagnitudeLOG',5) ;

if ~isempty(DATA.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG)
    Sincre = log10(S(1:end-1)./S(2:end)) ;
    
    III = find(Sincre >= DATA.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG) ;
    
    if ~isempty(III)
        nmodes= III(1) ;
        Lambda = Lambda(:,1:nmodes) ;
        S = S(1:nmodes) ;
        V = V(:,1:nmodes) ;
    end
    
    
end





% else
%
%     if  DATA.ndomainsSVD == 1
%         % dbstop('37')
%         [Lambda,S,~,ERRORsvd] = SVDtruncated(Xf,DATA.TOLsvdXf,0);
%     else
%         % dbstop('23')
%         DATA.PLOTsingVAL_Xf_All = 0 ;
%         DATA.TOLdomains = [] ;
%         [Lambda,S,~,ERRORsvd] = SVD_partROW(Xf,DATA.ndomainsSVD,DATA,[],size(Xf,1)) ;
%     end
% end

% if  DATA.SVD_BEFORE_zeroaverage ==1
%     % Basis vectors for Xf (W-orthogonal) --> inv(sqrt(W))*Lambda
%     % inv(sqrt(W))*(Lambda*S) plays the same role as Xf in the standard version
%     % Accordingly, we make
%     Xf = bsxfun(@times,Lambda',S)' ;
%     Winv = 1./sqrt(W) ;
%     Xf = bsxfun(@times,Xf,Winv) ;
%     % Then we calculate the zero.-average matrix
%     INTexac = Xf'*W ;
%     Xf = bsxfun(@minus,Xf',INTexac/V)' ;
%     Xf = bsxfun(@times,Xf,sqrt(W)) ;
%     % Anf finally, the SVD ...
%     disp('Additional SVD (for orthonormalizing )')
%      DATA.TOLsvdXf = [] ;
%     if  DATA.ndomainsSVD == 1
%         [Lambda,S,~,ERRORsvd] = SVDtruncated(Xf,DATA.TOLsvdXf,0);
%     else
%         % dbstop('23')
%         DATA.PLOTsingVAL_Xf_All = 0 ;
%         DATA.TOLdomains = [] ;
%         [Lambda,S,~,ERRORsvd] = SVD_partROW(Xf,DATA.ndomainsSVD,DATA,[],size(Xf,1)) ;
%     end
%
% end




if  DATA.PLOTsingVAL_Xf_All ==1
    figure(600)
    % dbstop('25')
    
    plot(ERRORsvd*100);
    xlabel('Mode Xfs ')
    ylabel('SVD error (%)')
    title('SVD error of Xf (for determining Jmin)')
    
    figure(601)
    % dbstop('25')
    %dbstop('53')
    plot(log10(ERRORsvd(find(ERRORsvd>0))));
    xlabel('Mode Xfs ')
    ylabel('log(SVD error) ')
    title('SVD error of Xf (for determining Jmin)')
end
%


% else
%     Lambda = Xf ; S =[] ;
% end


%[Lambda,S,~] = svd(Xf,0) ;
% -------------------------------
% S = diag(S) ;
% SingVsq = (S.*S) ;
% DENOM = sum(SingVsq) ;
%dbstop('168')
if DATA.includeSingularValuesinJ >0
    %dbstop('70')
    %     warning('borrar')
    %     S(1) = S(1)*1000 ;
    %      S(end) = 1000*S(1) ;
    factor = DATA.includeSingularValuesinJ  ;
    Lambda = bsxfun(@times,Lambda',(S).^factor)';
end



%dbstop('179')
p = size(Lambda,2) ;  % Number of modes
disp(['Number of selected modes (p) = ',num2str(p)])
% Matrix J appearing in the minimization problem  (step 5 in Box 5.1)
vol = sum(W);
if  DATA.ImposeVolumeConstraint == 1
    if DATA.bMIN_withONE ==1
        J = [Lambda';  sqrt(W')/vol] ;   %
        b = zeros(size(J,1),1) ;
        b(end) =   1;  %
    elseif  DATA.bMIN_withONE ==2
        J = [Lambda';  sqrt(W')/sqrt(vol)] ;   %
        b = zeros(size(J,1),1) ;
        b(end) =   sqrt(vol);  %
    else
        J = [Lambda';  sqrt(W')] ;   %
        b = zeros(size(J,1),1) ;
        b(end) =   V;  %
    end
else
    if  DATA.IncludeWeights_Jmin ==1
        J = Lambda' ;
        b = J*sqrt(W) ;
    else
        J = Lambda' ;
        b = J*(W) ;
    end
end
% For later purposes it is convenient to precompute the norm of J
Jnorm = sqrt(sum(J.*J,1)) ;

DATAOUT.SingV_Jmin = S ;
DATAOUT.V_Jmin = V ;
DATAOUT.INTexac = INTexac;
