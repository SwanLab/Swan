function [Q B]= RANDQB_adapt(C,R,mu,DATA)
%
% RANDQB_adapt iteratively construct an orthogonal basis matrix for the
% range of matrix C, such that  norm(C-Q*B,'fro') <= mu, where B = Q^T*C
%
% INPUTS:
% ------
% C --> M x N matrix
%
%  Optional inputs
%  ---------------
%  R is an estimation for an upper bound of rank(C)
% mu --> Tolerance (machine epsilon parameter)
%
% if nargin == 1, then R = ceil(DATA.rho_est*min(size(C))),
%    and mu =  min(size(C))*eps(norm(C,'fro'))
%
% DATA.NITER = Maximum number of iterations. Default = 20
%
% OUTPUT:
% ------
% Q --> Orthogonal basis matrix such that norm(C-Q*Q'*C,'fro') <= mu
% B = Q'*C
%
% This function is partially based on the randomized range finder algorithm
% proposed by Martisoon et al (2015) "A randomized blocked algorithm for efficiently computing
% rank-revealing factorizations of matrices".  The main novelty of our
% approach is the use of random matrices of varying size (the size of each random matrix is
%  estimated based on linear rank predictions).
% Written by Joaquín  A. Hernández, December 2016.
% UPC/CIMNE, Barcelona, Spain
% jhortega@cimne.upc.edu
%

%dbstop('33')
%dbstop('37')
if nargin == 0
    load('tmp.mat')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA = DefaultField(DATA,'RANKcubicINTERP',1) ;

[M,N] = size(C);
MNmin = min(size(C)) ;
%dbstop('48')
if nargin ==1
    R = ceil(0.05*min(size(C))); % Rank estimation (initial value)
    nC = norm(C,'fro') ; % Norm of the initial residual
    mu =   (max(size(C))*eps(nC))  ;  ; % MAchine precision parameter
    DATA = [] ;
elseif nargin ==2
    nC = norm(C,'fro') ; % Norm of the initial residual
    mu =   (max(size(C))*eps(nC))  ;  ; % MAchine precision parameter
    DATA = [] ;
elseif   nargin ==3
    DATA= [] ;
end

%dbstop('61')
if   isempty(mu) || mu==0
    nC = norm(C,'fro') ; % Norm of the initial residual
    mu =   (max(size(C))*eps(nC))  ;  ; % MAchine precision parameter
end



DATA = DefaultField(DATA,'NITER',10) ;
DATA = DefaultField(DATA,'Rest',[]) ;
DATA = DefaultField(DATA,'rho_est',0.05) ;

NITER = DATA.NITER ;

if R == 0
    R = ceil(DATA.rho_est*min(size(C))); % Rank estimation (initial value)
end

Omega = randn(N,R) ;
disp(['Testing rank =  ',num2str(R),' ...'])
Q = ORTH(C*Omega,mu) ; % Initial basis matrix
i = 1;
if ~isempty(Q)
    % If size(Q,2) < R, then Q is the desired basis matrix, and no more
    % operations are required. Otherwise, if
    Rhist =size(Q,2) ;
    if size(Q,2) == R
        % We have to incrementally augment Q
        %   dbstop('89')
        dRmin = ceil((MNmin-R)/NITER) ;  %  Minimum size for rank increment
        
        dR = dRmin ; % Size of the first rank increment
        B = Q'*C ;
        C = C - Q*B ; % Residual
        nC = norm(C,'fro') ; % Norm of the residual
        disp(['iter = ',num2str(0),'  nC =',num2str(nC),' R =',num2str(size(Q,2))])
        R = R + dR ;
        nChist =nC ;
        while nC>=mu
            disp(['Testing rank =  ',num2str(R),' ...'])
            nC_old = nC ;
            R_old = size(Q,2) ;
            Omega = randn(N,dR) ; %  Draw a N x dR random matrix
            %             if isempty(C) | isempty(Omega)
            %                 Omega ;
            %             end
            Qi = ORTH(C*Omega,mu); % Orthogonal basis matrix for residual
            if ~isempty(Q)
                Qi = orth(Qi - Q*(Q'*Qi)) ; % Reorthogonalization
            end
            Bi = Qi'*C ;
            C = C - Qi*Bi; % New residual
            Q = [Q Qi] ;  % Basis matrix is augmented with Qi
            B = [B; Bi];
            nC = norm(C,'fro') ; %
            disp(['iter = ',num2str(i),'  nC =',num2str(nC),' R =',num2str(size(Q,2))])
            R_new = size(Q,2);
            Rhist(i+1) =R_new ;
            nChist(i+1) = nC ;
            % Estimated rank (linear)
            %   if DATA.RANKcubicINTERP == 0 | i<3
            Rest = R_old+ (R_new-R_old)/(nC-nC_old)*(mu - nC_old) ;
            %  else
            %     RestLIN = R_old+ (R_new-R_old)/(nC-nC_old)*(mu - nC_old) ;
            %    Rest = interp1(nChist(end-3:end),Rhist(end-3:end),mu,'cubic','extrap') ;
            %end
            % Estimated increment
            dR = ceil(Rest-R_new) ;
            dR = max(dR,dRmin);
            R = R_new+ dR ;
            if R >=N
                %%%% The estimated rank is larger than the number of
                %%%% columns. It makes no sense to employ randomization to
                %%%% reveal the rank of C
                Q ='FULL' ;
                break
            end
            i = i + 1 ;
        end
        
        %   disp(['Upper bound for rank = ',size(Q,)])
    else
        disp(['Upper bound for rank = ',num2str(size(Q,2))])
        B = Q'*C ;
    end
else
    B= [] ;
end

%
% figure
% hold on
%
%     plot(Rhist,nChist,'r*')
%     xlabel('Rank')
%     ylabel('Norm')

