function [Lambda,INTexac,VSinv,A] = GetBasisMatrixECM(A,W,DATA)
% SVD of A*sqrt(W)
% --------------------------------------
%dbstop('5')
if nargin == 0
    load('tmp.mat')
    
end



% Exact integral
% ---------------
if ~iscell(A)
    INTexac = A'*W ; % High-fidelity integral of all snapshots
    A = bsxfun(@times,A,sqrt(W)) ;
else
    disp('Computing exact integrals, and multiplying by sqrt(W)')
    INTexac = cell(size(A)) ;
    for imat = 1:length(A)
        disp(['iblock = ',num2str(imat)])
        if isnumeric(A{imat})
            INTexac{imat} =A{imat}'*W ;
            A{imat} = bsxfun(@times,A{imat},sqrt(W)) ;
        else
            disp(['Retrieving from memory ...'])
            SSS = load(A{imat}) ;
            disp([' ... Done'])
            fff = fieldnames(SSS) ;
            Ai = SSS.(fff{1}) ; SSS = [] ;
            INTexac{imat} =Ai'*W ;
            Ai = bsxfun(@times,Ai,sqrt(W)) ;
            disp(['Saving again in  memory ... (multiplied by sqrt(W))'])
            save(A{imat},'Ai') ;
            disp([' ... Done'])
        end
    end
    INTexac = cell2mat(INTexac') ;
end


% if ~iscell(A)
%     A = bsxfun(@times,A,sqrt(W)) ;
% else
%     for ipart = 1:length(A)
%         if isnumeric(A{ipart})
%             A{ipart} = bsxfun(@times,A{ipart},sqrt(W)) ;
%         else
%             SSS = load(A{imat}) ;
%             fff = fieldnames(SSS) ;
%             %    Ai = SSS.(fff{1}) ; SSS = [] ;
%             INTexac{imat} =SSS.(fff{1})'*W ;
%         end
%     end
% end


% Orthogonal basis matrix Lambda for the column space of A
disp('SVD of Aw ')
DATA = DefaultField(DATA,'UsePartitionedRandomizedAlgorithm',1) ;

% if iscell(A)
%     DATA.UsePartitionedRandomizedAlgorithm = 1 ;
%
% end
DATA = DefaultField(DATA,'TOL_SVD_A',0) ;

if  ~iscell(A)
    if DATA.UsePartitionedRandomizedAlgorithm == 1
         DATAlocSVD.HIDE_OUTPUT = 0 ;
        tic
        [Lambda,S,V ] =      RSVDTrowblock({A},DATA.TOL_SVD_A,DATAlocSVD);
        toc
    else
        DATAsvd.RELATIVE_SVD =1 ;
        [Lambda,S,V] =    SVDT(A,DATA.TOL_SVD_A,DATAsvd ) ;
    end
else
    EPSILON_GLO = DATA.TOL_SVD_A*ones(size(A)) ;
    DATAlocSVD.HIDE_OUTPUT = 0 ;
    tic
    [Lambda,S,V ] =      RSVDTrowblock(A,DATA.TOL_SVD_A,DATAlocSVD );
    toc
    disp(['Number of modes =',num2str(length(S)),' (of ',num2str(size(V,2)),'columns'])
end
VSinv = bsxfun(@times,V',1./S)' ;


disp(['Number of left singular vectors  = ',num2str(size(Lambda,2))])

% % Enlarge the basis matris for Aw
a  = sqrt(W) - Lambda*(Lambda'*sqrt(W)) ;
TOL  = 1e-10 ; TOLmax = 0.9 ;
if  DATA.Method_Evaluation_Basis_Integrand == 2
    if norm(a) > TOLmax*norm(sqrt(W))
        error('Include in the list of integrand functions the constant functions (otherwise the problem will be  ill-posed)')
    end
    
else
    if norm(a) > TOL
        INCLUDE_ADDITIONAL_COLUMN = 1;
    else
        INCLUDE_ADDITIONAL_COLUMN = 0 ;
    end
    if    INCLUDE_ADDITIONAL_COLUMN ==1
        a = a/norm(a) ;
        Lambda = [a,Lambda] ;
    end
end



% IntApproxSVD = ErrorApproximateIntegral(A,Lambda,INTexac,W,DATA) ; 



% PLOT_SVDerror = 1;
% if  PLOT_SVDerror ==1
%     figure(356)
%     hold on
%     title('SVD error ')
%     SingVsq =  (S.*S) ;
%     SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
%     %normEf2 = sqrt(sum(SingVsq) - cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
%     normEf2 = sqrt(cumsum(SingVsq)) ;
%     normEf2 = sort(normEf2,'descend');
%     hold on
%     lS = log10(normEf2);
%
%     h = plot([0:length(S)-1],lS,'r');
%
%     xlabel('Number of Modes')
%     ylabel('log(error)')
%
%     %legend(h)
% end



%
% DATA = DefaultField(DATA,'TruncationCriterionSVDmodesJumpOrderMagnitudeLOG',5) ;
% if ~isempty(DATA.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG)
%     Sincre = log10(S(1:end-1)./S(2:end)) ;
%     III = find(Sincre >= DATA.TruncationCriterionSVDmodesJumpOrderMagnitudeLOG) ;
%     if ~isempty(III)
%         nmodes= III(1) ;
%         Lambda = Lambda(:,1:nmodes) ;
%         S = S(1:nmodes) ;
%         V = V(:,1:nmodes) ;
%     end
% end
%
% %%
%
%
