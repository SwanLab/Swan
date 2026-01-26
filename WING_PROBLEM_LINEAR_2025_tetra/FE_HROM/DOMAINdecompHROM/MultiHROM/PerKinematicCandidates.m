function WdefCAND = PerKinematicCandidates(faceDOFS,PhiDEF,Vrb,Mintf,TOLcosINTfluc)

if nargin == 0
    load('tmp1.mat')
end

% COSINE_ANGLES = [] ;
% idoms = find(I == cnINTF) ;
% [jdomGLO,jfaceGLO] =   ind2sub(size(cnINTF),idoms) ;
% 
% if length(jdomGLO) == 1
%     % End interface
%     % Extracting the fluctuating part of
%     fi = faceDOFS{jdomGLO}{jfaceGLO} ;
%     PhiDEF_f_i = PhiDEF{jdomGLO}(fi,:) ;
% %     WdefCAND =  PprojDEF_operator(Vrb{I},Mintf{I},PhiDEF_f_i);
%     
% else
%     %   \PhiFLUC{-}{I} \defeq    \PprojDEF{}{I} \Dintf{e^{-1}}{2} \PhiDEFf{e}{2},
    PhiFLUC = cell(1 ,2) ;
    for idomLOC = 1:2
%         INDEX_DOMAIN = jdomGLO(idomLOC) ;
%         INDEX_FACE =  idomLOC(idomLOC) ;
        fi = faceDOFS{idomLOC} ;
        PhiDEF_fi =  PhiDEF(fi,:) ;
        PhiFLUC{idomLOC} =   PprojDEF_operator(Vrb{idomLOC},Mintf{idomLOC},PhiDEF_fi);
    end
    % Compute the intersection between both subpsaces
    DATALOC.TOL = 1e-6; 
    TOLcosINTfluc = 1e-4 ; 
   % WdefCAND =   WSVDT([PhiFLUC{1},PhiFLUC{2}],Mintf{1},DATALOC) ; 
     [COSINE_ANGLES,WdefCAND]= PRANGLES(PhiFLUC{1},PhiFLUC{2},Mintf{1},1-TOLcosINTfluc,1e20) ;
%     
%     
% end