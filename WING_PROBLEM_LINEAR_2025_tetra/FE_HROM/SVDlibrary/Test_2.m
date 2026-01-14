clc
clear all

%%
RUN_CLASSIC = 0 ;
GENERATE_AGAIN = 1 ;
NAMEWS = 'DATAWS/data1.mat'
m = 2000; % Number of columns
n = m   ;  % Number of rows
r = 1000 ; % Estimated rank
e0 = 0 ;   % Relative tolerance
mu = 0 ;  % Relatine machine precision parameter
R =300;
DATA.NITER = 100 ;
DATA.TypeRankEstimate =0 ; % Exponential
%%%
SIZEM = m*n*8e-6 ;
disp(['SIZE = ',num2str(SIZEM),' MB'])

if GENERATE_AGAIN == 1
    disp('Generating matrix ...')
    % Orthogonal matrix of dimensions mxr
    U =randn(m,r) ;
    [U] = SVDT(U) ;
    % Orthogonal matrix of dimensions nxr
    V =randn(m,r) ;
    [V] = SVDT(V) ;
    % Singular values (exponential decay)
    % MAchine precision parameter
    mu = min(m,n)*eps ;
    % We set lambda so that the last SV is equal to mu
    lambda = -log(mu)/(r-1) ;
    nmodes = 1:r ;
    S =   exp(lambda*(1-nmodes))';
    figure(1)
    hold on
    xlabel('Mode')
    ylabel('Singular value')
    plot(nmodes,S)
    %%%%--------------------
    A = bsxfun(@times,V',S);
    A = U*A ; clear U; clear V
    disp(' ... End')
    disp('Saving matrix ...')
    save(NAMEWS,'A')
    disp('End matrix')
else
    disp('Loading matrix ...')
    load(NAMEWS,'A')
    disp('End matrix')
end


%%%%%%%%%%%%%%%%%%
% Standard SVD
%%%%%%%%%%%%%%%%%%%
if RUN_CLASSIC == 1
    disp('Standard SVD...')
    DATA.RELATIVE_SVD = 1;
    tic
    [~,S,~,eSVDclassic]  =SVDT(A,e0,DATA) ;
    disp('...End')
    TIME = toc ;
    LEGG = 'CLASS.SVD'
    eSVD = eSVDclassic ;
    LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVD),' time =',num2str(TIME),' s']
    figure(2)
    hold     on
    h =  plot(log10(S),'r')
    legend(h,LEGG)
end

%%%%%%%%%%%%%%%%%
% Randomized SVD
%%%%%%%%%%%%%%%%%
a = norm(A,'fro') ;
mu = mu*a ;
tic
e0 = e0*a ;
DATA.RELATIVE_SVD =0 ;
disp('Random. SVD...')
[~,S,~,eSVDrandom]  =RSVDT(A,e0,mu,R,DATA) ;
disp('...End')
TIME = toc
eSVD = eSVDrandom ;
LEGG = ['RSVDT'] ;
LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVD),' time =',num2str(TIME),' s']
figure(2)
hold     on
h =  plot(log10(S),'b--')
legend(h,LEGG)

%%%%%%%%%%%%%%%%%%%%%%%%%
legend off
legend show

