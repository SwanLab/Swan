function [DOFm,DOFs] = CoarseDOFS_DEIMbased...
    (D,DATAIN,Irb,nDOFsFACEall)

if nargin == 0
    load('tmp1.mat')
end

nfaces = length(nDOFsFACEall) ;
IndexesFaces = cell(nfaces/2,nfaces/2) ;
iini = 1;
for iface= 1:nfaces
    ifin = iini + nDOFsFACEall(iface) -1;
    IndexesFaces{iface} = (iini:ifin)' ;
    iini = ifin +1 ;
end
DOFsCOARSE = cell2mat(IndexesFaces) ;  % First column --> Faces 1 and 2. Second column --> Faces 3 and 4

IrbPAIRS = zeros(length(Irb)/2,2) ;
IrbPAIRS(:,1) = sort(intersect(DOFsCOARSE(:,1),Irb)) ;
IrbPAIRS(:,2) = sort(intersect(DOFsCOARSE(:,2),Irb)) ;


% STEP 1
% ------
% Initial DATA: $\D$ and $I_{rb}$. Necessary condition: $\D(I_{rb},:)$ must be full rank.
Dirb = D(Irb,:) ;
if rank(Dirb)  ~= size(Dirb,1)
    error('The initial set of rigid body DOFs is conducive to ill-posed problems')
end

% STEP 2
% ---------------
% Select a set of columns  $J$ such that  $\D(I_{rb},J)$ is invertible.
% This can be done by applying the DEIM
%   \begin{equation}
%    J = DEIM(\D(I_{rb},:))
%   \end{equation}
J = DEIMfun(Dirb') ;

% STEP 3
% ----------------------------
% Let  $L$ be the complementary set of $J$, and let us set $k=0$,  $\A = \D(:,J)$    and $I = I_{rb}$.
L = 1:size(D,2) ; L(J)= [] ;
A = D(:,J) ; I = Irb ;
IcandidatePAIRS = DOFsCOARSE;
IcandidatePAIRS(IrbPAIRS(:,1),:) = [] ;

for kLOC=1:(length(L)/2)
    %  1) Compute the residual
    % \begin{equation}
    %  \r = \D(:,L(k)) - \A \A(I,:)^{{-1}}  \D(:,L(k))
    % \end{equation}
    k = 2*kLOC-1:2*kLOC ;
    r = D(:,L(k)) - A*(A(I,:)\D(I,L(k)));
    r = sqrt(sum(r.^2,2)) ;
    % Notice that we are trying to write the next column of $\D$ as a linear combination of the already selected columns ($\A$) ---via an \emph{oblique} projection.
    % \item STEP 2) Criterion for choosing the best next pair of rows.
    % We partition the residual as  $\r = \coldos{\r_{-}}{\r_{+}}$, and calculate
    % \begin{equation}
    %  a_i = \absval{\r(\DOFsCOARSE{-}{}(i))} \absval{\r{(\DOFsCOARSE{+}{}(i))}},   \hspace{0.5cm}  i = 1,2 \ldots \nDOFCdom{}{}/2
    % \end{equation}
    a = sqrt(r(DOFsCOARSE(:,1))).^2  +  (r(DOFsCOARSE(:,2)).^2) ;
    
    % The next row to incorporate into the set of master DOFs is that maximizing $a_i$.
    % \begin{equation}
    %  i = \argMAX{a_i}
    % \end{equation}
    [aMAX,iCANDIDATE]= max(a(IcandidatePAIRS(:,1))) ;
    %i = IcandidatePAIRS(iCANDIDATE,1) ;
    
    % Set  $I \leftarrow I \cap \DOFsCOARSE{+}{i} \cap \DOFsCOARSE{-}{i} $
    % and $\A \leftarrow [\A,\D(:,L(k))]$,  $k \leftarrow k +1$. Go to STEP 1
    Inew = IcandidatePAIRS(iCANDIDATE,:) ;
    
    I = [I Inew] ;
    A = [A,D(:,L(k))];
    
    IcandidatePAIRS(iCANDIDATE,:) =[] ;
    
    
end

DOFm = sort(I) ;
DOFs = 1:size(D,1) ;
DOFs(DOFm) =[];


DATAIN.KINEMATIC_CONSTRAINTS_MODES = ...
    DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'INCLUDE_SLAVE_DOFs',0)  ;
if DATAIN.KINEMATIC_CONSTRAINTS_MODES.INCLUDE_SLAVE_DOFs == 0
    DOFs=[];
end

SSS =svd(D(DOFm,:)) ;

if rank(D(DOFm,:)) ~= length(DOFm)
    error('Ill-posed problem....')
end

%
% % STEP 6. Master and Slave DOFs
% % Select linearly independent columns
% % --------------------------------------
% % First we have to select the candidates for being master DOFs
% % The first criterion is the angle formed by mode  I+  with mode I-
% % Only those angles above a certain threshold are selected
% % Dcomp is a nDOFT x nDOFM matrix. The goal is to select a set of rows  DOFm   so that Dcomp(DOFm,:) is invertible
% % However, the selection should be done observing the conditions of
% % periodicity. Indeed, the candidate modes are ordered in pairs.  We  have
% % to ascertain how linearly correlated are these pair of modes. In the
% % following function, we calculate the angle formed by each pair of modes.
% angCOV = CorrelationAnglesPairFaces(Dcomp)  ;
% % New we select only those angles that are above a certain threshold
% DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'ANGLE_THRESHOLD_COVARIANCE',10) ;  ;
% angTOL = DATAIN.KINEMATIC_CONSTRAINTS_MODES.ANGLE_THRESHOLD_COVARIANCE ; % Tolerance for filtering out modes
%
% METHOD_SELECT = 1;
%
% if METHOD_SELECT == 1
%     % Simplified method
%      nmodesP = size(BasisUdom_f,2)/2 ;  % Number of required modes (for one side)
%
%      DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField...
%          (DATAIN.KINEMATIC_CONSTRAINTS_MODES,'RIGID_BODY_MODES_TO_INCLUDE',{[1 2],[1,2]})  ;
%
%      if ~isempty(DATAIN.KINEMATIC_CONSTRAINTS_MODES.RIGID_BODY_MODES_TO_INCLUDE)
%          iRB = DATAIN.KINEMATIC_CONSTRAINTS_MODES.RIGID_BODY_MODES_TO_INCLUDE ;
%          % Indices for rB side 3
%          iRB1 = intersect(iRB{1},1:nRB(1)) ;
%        %  iRB1 =iRB{1} ;
%                 iRB2 = intersect(iRB{2},1:nRB(1)) ;
%
%
%          iRB2 =  nDOFsFACEall(1)+iRB2 ;
%
%          IndMaster_m = [iRB1,iRB2]' ;
%
%          % Remaining entries
%          nREST = nmodesP-length(IndMaster_m) ;
%          if nREST > 0
%              indDEF = setdiff(1:length(angCOV),IndMaster_m) ;
%              angCOVdef = angCOV(indDEF) ;
%                [MMM,III]  = sort(angCOVdef,'descend') ;
%                 III_rest = III(1:nREST) ;
%                 IndMaster_mrest = indDEF(III_rest) ;
%                 IndMaster_m =[IndMaster_m; IndMaster_mrest(:)] ;
%          end
%
%      else
%          [MMM,III]  = sort(angCOV,'descend') ;
%
%          IndMaster_m = III(1:nmodesP) ;
%      end
%
%
%
%
%      IndMaster_p = IndMaster_m +  size(Dcomp,1)/2 ;
%      DOFm = sort([IndMaster_m;IndMaster_p]) ;
%
%      %    DOFm = IndicesRB ;
%
% else
%
%     TTT=  find(MMM > angTOL) ;
%     % Candidates for being master DOFs (negative  sides, i.e, 1-2)
%     CandLOC_m = III(TTT);
%
%     if length(CandLOC_m)*2 < size(BasisUdom_f,2)
%         error('SET ANGLE_THRESHOLD_COVARIANCE to a LOWER VALUE, or perhaps reduce the number of reaction modes')
%     end
%
%     % Now we have to select among CandLOC_m and
%     CandLOC_p = CandLOC_m +  size(Dcomp,1)/2 ;
%
%
%     % most linearly independent rows of Dcomp
%     INDICES_SELECTED = [1] ;
%     nmodesP = size(BasisUdom_f,2)/2 ;  % Number of required modes (for one side)
%     DnT = Dcomp(CandLOC_m,:)' ;
%     for icand = 2:nmodesP
%         IndRows = setdiff(1:size(DnT,2),INDICES_SELECTED) ;
%         OrthCompl = DnT(:,IndRows) - DnT(:,INDICES_SELECTED)*(DnT(:,INDICES_SELECTED)\DnT(:, IndRows) ) ;
%         normORTH = sqrt(sum(OrthCompl,1).^2) ;
%         [~,imin] = min(normORTH) ;
%         INDICES_SELECTED(icand) = IndRows(imin) ;
%     end
%
%     IndMaster_m = CandLOC_m(INDICES_SELECTED) ;
%     IndMaster_p = CandLOC_p(INDICES_SELECTED) ;
%
%     DOFm = sort([IndMaster_m;IndMaster_p]) ;
%
%
%
% end
%
%     DOFs = 1:size(Dcomp,1) ;
%     DOFs(DOFm) = [] ;
%
%   %  Dcomp_M  =Dcomp(DOFm,:) ;
%