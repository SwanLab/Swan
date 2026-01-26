clc
clear all

%%
TYPE_SV =  'EXP' ; 'POLY'  ;
RUN_CLASSIC = 0 ;
GENERATE_AGAIN = 1 ;
NAMEWS = 'DATAWS/data1.mat'
m = 3000; % Number of columns
n = m   ;  % Number of rows
r = 100 ; % Estimated rank
e0 = 0 ;   % Relative tolerance
mu = 0 ;  % Relatine machine precision parameter
R =50;
DATA.TypeRankEstimate =1 ; % Exponential
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
    switch  TYPE_SV
        case  'EXP'
            % We set lambda so that the last SV is equal to mu
            lambda = -log(mu)/(r-1) ;
            nmodes = 1:r ;
            S =   exp(lambda*(1-nmodes))';
        case 'POLY'
            p = 5 ;
            a = (1-mu)/((r-1)^p) ;
            nmodes = 1:r ;
            S = 1- a*(nmodes-1).^p ;
            S = S' ;
    end
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
    
    
    
    figure(3)
    hold on
    
    SingVsq =  (S.*S) ;
    SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
    %normEf2 = sqrt(sum(SingVsq) - cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
    normEf2 = sqrt(cumsum(SingVsq)) ;
    normEf2 = sort(normEf2,'descend');
    hold on
    lS = log10(normEf2);
    
    plot([0:length(S)-1],lS);
    
    xlabel('Modes')
    ylabel('log(error)')
    
    
    
       figure(4)
    hold on
    
 
    plot([0:length(S)-1],normEf2);
    
    xlabel('Modes')
    ylabel('error')
    
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
    DATA.COMPUTE_V = 0 ;
    tic
    [~,S,~,eSVDclassic]  =SVDT(A,e0,DATA) ;
    disp('...End')
    TIME = toc ;
    LEGG = 'CLASS.SVD';
    eSVD = eSVDclassic ;
    LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVD),' time =',num2str(TIME),' s']
    figure(2)
    hold     on
    h =  plot(log10(S),'r');
    legend(h,LEGG);
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
LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVD),' time =',num2str(TIME),' s'];
figure(2)
hold     on
h =  plot(log10(S),'b--');
legend(h,LEGG)
xlabel('Modes')
ylabel('Log(S)')

%%%%%%%%%%%%%%%%%%%%%%%%%
legend off
legend show

