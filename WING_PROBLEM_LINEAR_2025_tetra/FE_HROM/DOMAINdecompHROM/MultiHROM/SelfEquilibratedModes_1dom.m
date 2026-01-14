function  [PsiDEFf,MdomFFinv,AsnapREACse,MdomFFinv_chol] =   SelfEquilibratedModes_1dom(faceDOFS,PsiRBf,MdomFF,AsnapREAC,DATALOC)
% COMPUTATION OF SELF-EQUILIBRATED MODES 
% 


if nargin == 0
    load('tmp.mat')
elseif nargin == 4
    DATALOC = [] ; 
end
  
%   MdomFF = speye(size(MdomFF)) ; 
% disp('borrar, temporal')


[AsnapREACse,MdomFFinv] =  HprojSEf_operator(PsiRBf,MdomFF,AsnapREAC(faceDOFS,:)) ;
% Weighted SVD

  %DATALOC.TOL = 1e-12 ;
  
  
  
  DATALOC = DefaultField(DATALOC,'TOL',1e-6);
  
  
    [PsiDEFf,S,V,MdomFFinv_chol ]= WSVDT(AsnapREACse,MdomFFinv,DATALOC) ;
    
    
    
% if DATAinputLOC.SVD_PER_TRAINING_TEST_SELFEQ == 1
%     for itest = 1:ntests
%         AsnapREACse_test{itest} = AsnapREACse(:,itest:ntests:end) ;
%         nnn= norm(AsnapREACse_test{itest},'fro') ;
%         AsnapREACse_test{itest}  = AsnapREACse_test{itest}/nnn ;
%         
%         %     [~,SS,~]  = SVDT(AsnapREACse_test{itest}) ;
%         %     SS = SS/SS(1) ;
%     end
%     if length(TOLse) ==1
%         epsilon = ones(1,ntests)*TOLse ;
%     else
%         epsilon = TOLse ;
%     end
%     DATALOC = [] ;
%     [U,S,V,ETIME,eSVD,RankMatrix,DATAOUT] = RSVDqp(AsnapREACse_test,epsilon,DATALOC) ;
%     [PsiDEFf,S,V ]= WSVDT(U,cell2mat(MdomFFinv)) ;
%     
% else
  
%end






