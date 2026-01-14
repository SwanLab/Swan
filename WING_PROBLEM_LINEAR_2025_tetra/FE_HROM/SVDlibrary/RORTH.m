function [Q B nC c]= RORTH(C,mu,R,DATA)
% RORTH iteratively constructs an orthogonal basis matrix Q for the
% range of matrix C such that  norm(C-Q*B,'fro') <= mu, where B = Q^T*C,
% and Q^T*Q = I.
%
% Mandatory arguments:
%
% C --> M x N matrix
%
% Optional arguments
%
% R: estimation for  rank(C). Default value: R = ceil(DATA.rho_est*min(size(C)))
% mu:  Tolerance.  Default value: mu =   (max(size(C))*eps(nC))  where   nC = norm(C,'fro')
% DATA.NITER = Maximum number of iterations. Default = 10
% DATA.rho_est = Parameter defining the rank estimation R. Default value =
% 0.05
%
% OUTPUT:
% ------
% Q --> Orthogonal basis matrix such that norm(C-Q*Q'*C,'fro') <= mu
% B = Q'*C
%
% REMARK: If the estimated rank turns out to be greater than the number of
% columns of C, then randomization makes no sense (it is a full rank
% matrix) and the algorithm returns Q = 'FULL' and  B = [] ;
%
% This function is partially based on the randomized range finder algorithm
% proposed by Martisoon et al (2015) "A randomized blocked algorithm for efficiently computing
% rank-revealing factorizations of matrices".  The main novelty of our
% approach is the use of random matrices of varying size (the size of each random matrix is
%  estimated based on linear rank predictions).
%
% Written by Joaquín  A. Hernández, March 2017.
% UPC/CIMNE, Barcelona, Spain
% jhortega@cimne.upc.edu
%dbstop('37')
if nargin == 0
    %     C = rand(1000,50) ;
    %     C =[C,C.^2,1000*C.^3] ;
    %     C = [C,6*C] ;
    %     C = [C;C] ;  C = [C C.^2] ;
    %       R = 0; % Rank estimation (initial value)
    %     nC = norm(C,'fro') ; % Norm of the initial residual
    %     mu =   (max(size(C))*eps(nC))  ;  ; % MAchine precision parameter
    %     DATA = [] ;
    load('tmp.mat')
    
    %  mu = 1e10 ;
end
[M,N] = size(C); MNmin = min(size(C)) ;
% ------------------------------------
% Default value for optional arguments
% ------------------------------------
if nargin ==1
    R = 0; % Rank estimation (initial value)
    nC = norm(C,'fro') ; % Norm of the initial residual
    mu =   (max(size(C))*eps(nC))  ;  ; % MAchine precision parameter
    DATA = [] ;
elseif nargin ==2
    nC = norm(C,'fro') ; % Norm of the initial residual
    mu =   (max(size(C))*eps(nC))  ;  ; % Machine precision parameter
    DATA = [] ;
elseif   nargin ==3
    DATA= [] ;
end
c = norm(C,'fro') ; nC = c ;  % Norm of the initial residual
if  isempty(mu) || mu==0
    mu =   (max(size(C))*eps(nC))  ;  ; % aAchine precision parameter
end
dRmax = ceil(0.25*min(size(C))) ;
dRmin = min(1,min(size(C))) ;
dRmin = max(dRmin,ceil(0.005*min(size(C)))) ;
DATA = DefaultField(DATA,'dRmax',dRmax) ; %  Maximum size   (relative)
DATA = DefaultField(DATA,'dRmin',dRmin) ; %  Maximum size   (relative)
DATA = DefaultField(DATA,'R',ceil(0.05*min(size(C)))) ; %  Maximum size   (relative)
DATA = DefaultField(DATA,'TypeRankEstimate',1) ;  % = 1, exponential
DATA = DefaultField(DATA,'PLOT_POINTS',0) ;  % = 1, exponential

DATA = DefaultField(DATA,'HIDE_OUTPUT',1) ;

%NITER = DATA.NITER ;
if R == 0
    R =DATA.R;
end
%dRmin = ceil((MNmin-R)/NITER) ;  %  Minimum size for rank increment
% -------------------
% End Default Values
% -------------------
Q = [] ; B = [] ; dR = R ; nC = c; i = 1; nC_old = c ; R_old = 0 ;
%dbstop('80')
DATA.COMPUTE_V = 0 ; DATA.RELATIVE_SVD = 0 ;
%dbstop('90')
%EXIT  = 0 ;
PLOT_POINTS = DATA.PLOT_POINTS;
PLOTFIGURE =3 ; 
if PLOT_POINTS == 1
    figure(PLOTFIGURE)
    hold on
    plot(R_old,log10(nC_old),'x')
    
%   HHN =   plot([0,100],[log10(mu),log10(mu)],'k')   ;
%   legend(HHN,['\mu'])
    
end
%dbstop('102')
while nC>=mu  %& EXIT ==0
    
    
    
    Omega = randn(N,dR) ; %  Draw a N x dR random matrix
    nOmega = sqrt(prod(size(Omega))) ;
    %DATA.SVD.ToleranceCutOff = 1 ;
    %   dbstop('90')
    %   [Qi] = SVDT(C*Omega/nOmega,mu,DATA); % Orthogonal basis matrix for residual
    factorRED = 10 ;
    DATA.SVD.MaxSize = max(size(C))/factorRED;
    %  [Qi] = SVDT(C*Omega/nOmega,0,DATA); % Orthogonal basis matrix for residual
    [Qi] = SVDT(C*Omega/nOmega,0,DATA); % Orthogonal basis matrix for residual
    %    if size(Qi,2) < dR ; EXIT = 1 ; end
    
    if isempty(Qi); break  ; end ;
    if ~isempty(Q)
        DATA.SVD.MaxSize = max(size(Qi));
        [Qi ]= SVDT(Qi - Q*(Q'*Qi),0,DATA) ; % Reorthogonalization
    end
    Bi = Qi'*C ; C = C - Qi*Bi; % Computing residual
    Q = [Q Qi] ;   B = [B; Bi]; % Basis matrix is augmented with Qi
    nC = norm(C,'fro') ; % Norm of the residual
    if DATA.HIDE_OUTPUT == 0
    disp(['iter = ',num2str(i),'  nC =',num2str(nC),'  dR =',num2str(dR),' R =',num2str(size(Q,2))]) 
    end
    
    %   if  i == 1  & DATA.TypeRankEstimate ==0
    %        dR = dRmin ; Rest = R + dR ;
    %   dR = ceil(0.1*R) ; Rest = R + dR ;
    %   else
    R_new = size(Q,2);
    if DATA.TypeRankEstimate == 0
        Rest = R_old+ (R_new-R_old)/(nC-nC_old)*(mu - nC_old) ;  % Estimated rank (linear)
    elseif DATA.TypeRankEstimate == 1
        Rest = R_old + (R_new -R_old)/(log(nC) - log(nC_old))*(log(mu) - log(nC_old)) ;
    else
        Rest = R_old + DATA.dRmin ; 
    end
    dR = ceil(Rest-R_new) ;
    dR = min(dR,DATA.dRmax);
    dR = max(dR,DATA.dRmin) ;
    Rest = R_new+ dR ;
    if Rest >=N ; Q ='FULL' ; B = [] ; ;   break ;  end
    %  end
    i = i + 1 ;
    
    R_old = R ; nC_old = nC ; R = Rest ;
    if PLOT_POINTS == 1
        figure(PLOTFIGURE)
        hold on
        plot(R_old,log10(nC_old),'ko')
    end
end


if PLOT_POINTS == 1
    figure(PLOTFIGURE)
    hold on
    
    aaa = axis ; 
     
  HHN =   plot([aaa(1),aaa(2)],[log10(mu),log10(mu)],'k')   ;
  legend(HHN,['\mu'])
     
end
