function  [PhiDEF,PsiDEFf,MdomCHOL] =   DeformationalModes_1dom(faceDOFSall,PsiRBf,MdomFF,BasisU,PsiDEFf,MdomFFinv,PhiRB,Mdom)  ;
% Computation of deformational modes
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% See explanation how it works in 
% 
if nargin == 0
    load('tmp.mat')    
end
%
%
% AsnapDISP_e = [] ;
% for eLOC = 1:length(domsSELECT)
%     e = domsSELECT(eLOC) ;
%     % Deformational  component of the snapshot matrix (projection operator)
%
%     if DATAinputLOC.NormalizeSnapshotsDEF == 1
%         normA = sqrt(sum(AsnapDISP{e}.^2,1))  ;
%         AsnapDISP_loc = bsxfun(@times,AsnapDISP{e}',1./normA')' ;
%     end
%     AsnapDISP_e = [AsnapDISP_e,  AsnapDISP_loc ] ;
%
%
% end
AsnapDISPdef =  SprojDEF_operator(PhiRB,Mdom,BasisU) ;
%
% ntests = size(AsnapDISP{1},2) ;
% AsnapDISPdef_test  = cell(1,ntests) ;
% if DATAinputLOC.SVD_PER_TRAINING_TEST_SELFEQ == 1
%     for itest = 1:ntests
%         AsnapDISPdef_test{itest} = AsnapDISPdef(:,itest:ntests:end) ;
%         nnn= norm(AsnapDISPdef_test{itest},'fro') ;
%         AsnapDISPdef_test{itest}  = AsnapDISPdef_test{itest}/nnn ;
%
%              [~,SS,~]  = SVDT(AsnapDISPdef_test{itest}) ;
%             SS = SS/SS(1) ;
%     end
% %     if length(TOLse) ==1
% %         epsilon = ones(1,ntests)*TOLse ;
% %     else
% %         epsilon = TOLse ;
% %     end
% epsilon= zeros(1,tests) ;
% DATALOC = [] ;
% [U,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(AsnapREACse_test,epsilon,DATALOC) ;
% [PsiDEFf,S,V ]= WSVDT(U,cell2mat(MdomFFinv)) ;
%
% else
%     DATALOC.TOL = TOLse ;
%     [PsiDEFf,S,V ]= WSVDT(AsnapREACse,cell2mat(MdomFFinv),DATALOC) ;
% end
f = faceDOFSall ;
%
% NEWMeTHOD = 1;
%
% if NEWMeTHOD ==1
% STEP 1
[y,S,V ]= SVDT(AsnapDISPdef(f,:),0) ; % Basis matrix for    deformational displacements (boundary)

coeffBASIS_y = bsxfun(@times,V',S)' ;
% coeffBASIS_y = bsxfun(@times,V',1./S)' ;  % WRONG !!!!  Amendment made on
% April 4th 2023...If S is not truncated, large values may appear. It is better to leave the "original" version, 
% appear in giuliodori2023multiscale.pdf 

% STEP 2
BASISFOUND =0  ;
while BASISFOUND == 0
    COEFF = y'*PsiDEFf ; % PROJECTION REACTION MODES ONTO SPAN DEFORMATIONAL MODES
    % STEP 3. Coefficients projection
    [UU,SS,VV] = SVDT(COEFF);  % SVD coefficients
    TOL = 1e-10 ;
    if find(SS/SS(1) <1e-10)
        error('This problem is not amenable to reduction. There are buble modes  ')
    end
    % STEP 4
    [PhiDEF,~,~,MdomCHOL] = WSVDT(AsnapDISPdef*(coeffBASIS_y*UU),Mdom) ;
    
    if size(PhiDEF,2) < size(PsiDEFf,2)
        %                 % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DOMAINdecompositionGEN/TruncationAnalysis.mlx
        %                 % Recomputing PsiDEFf{e}
        %                 [SS,PsiNEW] = PRANGLES(PsiDEFf{e},PhiDEF{e}(f,:),[],0,1e20) ;
        %                  [PsiDEFf{e},S,V ]= WSVDT(PsiNEW,cell2mat(MdomFFinv{e})) ;
        PsiDEFf = PsiDEFf(:,1:size(PhiDEF,2)) ;
    else
        BASISFOUND = 1;
    end
    
end

% else
%
%
%     [BasisUnew]= ReactionDispModesEqual2022(AsnapDISPdef,PsiDEFf,f,DATAROM);
%     [PhiDEF,~,~] = WSVDT(BasisUnew,Mdom) ;
%
% end

% 
% 
% DATAinputLOC = DefaultField(DATAinputLOC,'PLOT_EVOLUTION_def_MODES_FOR_EACH_TEST',0) ;
% 
% if DATAinputLOC.PLOT_EVOLUTION_def_MODES_FOR_EACH_TEST
%     
%     ifigureBASE = 799;
%     
%     
%     ntests = size(AsnapDISP{1},2) ;
%     
%     ndofs_red = size(PhiDEF,2) ;
%     nslices = length(AsnapDISP) ;
%     for itests = 1:ntests
%         SNAPtests_proyect = zeros(ndofs_red,nslices);
%         for eLOC = 1:length(AsnapDISP)
%             SNAPtests_proyect(:,eLOC) = PhiDEF'*Mdom*AsnapDISP{eLOC}(:,itests) ;
%         end
%         
%         figure(ifigureBASE+itests)
%         hold on
%         title(['modal contribution (def.) to TEST = ',num2str(itests)]);
%         xlabel('Slice')
%         ylabel('Modal contribution')
%         h = zeros(ndofs_red,1) ;
%         LLL = cell(ndofs_red,1) ;
%         
%         for iplot = 1:length(h)
%             h(iplot) = plot(SNAPtests_proyect(iplot,:)) ;
%             LLL{iplot} = ['Mode =',num2str(iplot)] ;
%         end
%         legend(h,LLL)
%         
%         
%         
%     end
%     
%     
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %