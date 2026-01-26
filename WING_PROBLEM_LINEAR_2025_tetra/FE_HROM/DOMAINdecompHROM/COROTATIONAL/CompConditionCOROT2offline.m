function [YcmpD,Zfict] = CompConditionCOROT2offline(ndim,Hrb,Grb,PsiSE,PsiRES,RrbB,Hdef,PdownsDEFb,Trb,HrbBUB)
% REduced-order operators derived from the compatibility condition
% associated to rotations
% JAHO, 1-Dec-2024, Sunday. Balmes 185, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/04_UPDrotUNC.mlx
% Latex notation:: /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTexAPPV.tex
if nargin == 0
    load('tmp1.mat')
end


if ndim == 2
    TRind = 1:2  ;
    ROTind = 3;
    identLOC = 1 ;
elseif ndim == 3
    TRind = 1:3  ;
    ROTind = 4:6;
    identLOC = eye(3) ;
else
    error('Option not implemented')
end

% \Hrb  =    \PsiRES^T \PhiDEFb = \coldos{\HrbTR}{\HrbROT} = \coldos{\PsiREStr^T\PhiDEFb}{\PsiRESrot^T \PhiDEFb}
%\GrbTR \defeq  \PsiRES^T \RrbTRb =  \coldos{\GrbTRtr}{\GrbROTtr} =      \coldos{\PsiREStr^T \RrbTRb}{\PsiRESrot^T \RrbTRb}
HrbTR = Hrb(TRind,:) ;
HrbROT = Hrb(ROTind,:) ;
GrbTRtr = Grb(TRind,TRind) ;
GrbROTtr = Grb(ROTind,TRind) ;

if ndim == 2
    % \vMIXe{e} =  (\GrbROTtr    \GrbTRtr^{-1})^{(e)}
    vMIX  = GrbROTtr*inv(GrbTRtr) ;
    % \jDEFe{e} \defeq   \Par{ \Par{ \vMIX  \HrbTR  - \HrbROT }  \Hdef^{-1}   }^{(e)}
    jDEF = (vMIX*HrbTR-HrbROT)*inv(Hdef) ;
    % \jRBe{e}  \defeq    \rowdos{-\vMIXe{e}}{\ident}
    jRB = [-vMIX,identLOC] ;
    % \jFICTe{e} \defeq  \rowdos{\jRBe{e}}{\jDEFe{e}}
    jFICT = [jRB,jDEF] ;
    % -------------------------------
    
    %  \cDEFrE{e}{i} \defeq  \colcuatro{\PsiSEe{e}(1:2:\MdofsF,i)^T \RrbROTbE{e}(1:2:\MdofsF) }
    %  {\PsiSEe{e}(2:2:\MdofsF,i)^T \RrbROTbE{e}(1:2:\MdofsF)}
    %  {\PsiSEe{e}(1:2:\MdofsF,i)^T \RrbROTbE{e}(2:2:\MdofsF)}
    %  {\PsiSEe{e}(2:2:\MdofsF,i)^T \RrbROTbE{e}(2:2:\MdofsF)},  \hspace{1cm}  i = 1,2 \ldots \pDEF
    RrbROTb = RrbB(:,ROTind) ;
    cDEFr = zeros(4,size(PsiSE,2)) ;
    for idef = 1:size(PsiSE,2)
        cDEFr(1,idef) = PsiSE(1:2:end,idef)'*RrbROTb(1:2:end) ;
        cDEFr(2,idef) = PsiSE(2:2:end,idef)'*RrbROTb(1:2:end) ;
        cDEFr(3,idef) = PsiSE(1:2:end,idef)'*RrbROTb(2:2:end) ;
        cDEFr(4,idef) = PsiSE(2:2:end,idef)'*RrbROTb(2:2:end) ;
    end
    
    %
    % \cRBtrE{e}{i} \defeq  \colcuatro{\PsiREStrE{e}(1:2:\MdofsF,i)^T \RrbROTbE{e}(1:2:\MdofsF) }
    % {\PsiREStrE{e}(2:2:\MdofsF,i)^T \RrbROTbE{e}(1:2:\MdofsF)}
    % {\PsiREStrE{e}(1:2:\MdofsF,i)^T \RrbROTbE{e}(2:2:\MdofsF)}
    % {\PsiREStrE{e}(2:2:\MdofsF,i)^T \RrbROTbE{e}(2:2:\MdofsF)},
    
    %\cRBrE{e} =  \rowtres{\cRBtrE{e}{1}}{\cRBtrE{e}{2}}{\cRBrrE{e}}   
    
    cRBr = zeros(4,size(PsiRES,2)) ;
    for idef = 1:size(PsiRES,2)
        cRBr(1,idef) = PsiRES(1:2:end,idef)'*RrbROTb(1:2:end) ;
        cRBr(2,idef) = PsiRES(2:2:end,idef)'*RrbROTb(1:2:end) ;
        cRBr(3,idef) = PsiRES(1:2:end,idef)'*RrbROTb(2:2:end) ;
        cRBr(4,idef) = PsiRES(2:2:end,idef)'*RrbROTb(2:2:end) ;
    end
    
    %  \cRe{e} \defeq  \rowdos{\cRBrE{e}{}}{\cDEFrE{e}{}}
    cR = [cRBr,cDEFr] ;   
    % \ZfictE{e}  =  {\jFICTe{e}}  (\cRe{e})^T
%     Zrb = jRB*cRBr' ; 
%     Zdef = jDEF*cDEFr' ; 
      Zfict = jFICT*cR' ; 
    
else
    error('3D option not implemented yet')
end

TrbROT = Trb(ROTind,:) ; 
TrbTR = Trb(TRind,:) ; 
HrbBUBrot = HrbBUB(ROTind,:) ; 
HrbBUBtr = HrbBUB(TRind,:) ; 

% Residual for the compatibility equation 

% 
%   \YcmpDe{e}  \defeq  \rowdos{\JrbROTmixE{e}}{\HrbBUBrotMIXe{e}}


% \HrbBUBrotMIX =  \HrbBUBrot -\GrbROTtr \GrbTRtr^{-1}\HrbBUBtr
 HrbBUBrotMIX = HrbBUBrot-GrbROTtr*(GrbTRtr\HrbBUBtr) ; 
% \JrbTR \defeq  \HrbTR \PdownsDEFb    - \TrbTR
JrbTR = HrbTR*PdownsDEFb-TrbTR ; 
% \JrbROT \defeq \HrbROT \PdownsDEFb    - \TrbROT
JrbROT = HrbROT*PdownsDEFb-TrbROT ; 
%   \JrbROTmix = \JrbROT  -\GrbROTtr \GrbTRtr^{-1}\JrbTR
JrbROTmix = JrbROT - GrbROTtr*(GrbTRtr\JrbTR) ; 
% \YcmpDe{e}  \defeq  \rowdos{\JrbROTmixE{e}}{\HrbBUBrotMIXe{e}}
YcmpD = [JrbROTmix,HrbBUBrotMIX] ; 
