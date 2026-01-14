clc
clear all
format long g

%
DATA.PLOT_POINTS = 1;  % Plot points arising from RANDOMIZATION
RUN_CLASSIC = 0;  % Run classical SVD
GENERATE_AGAIN = 0 ; % Generate again input data 
e0 = 0.001 ;   % Relative tolerance --- To truncate the SVD
mu = 0.00 ;  % Relatine machine precision parameter ---threshold for truncating the Q-search algorithm
NAMEDATA = 'DATA101' ; % General input data (to generate the function)
R =5;  % Initial rank estimate 
DATA.TypeRankEstimate =1; % Exponential == 1
DATA.dRmin = 1 ; % Minimum size of dR 
%%
addpath('AUXFUNCTIONS')
addpath('DATAINPUT')
eval(NAMEDATA)

if GENERATE_AGAIN == 1
    [A,alpha,beta,epsilon,p,q] =GenerateMatrixGlo(N,q0,M,p0,DATA,epsilon,nrepROWS,nrepCOLS,...
        COMPUTE_TRANSPOSE,DATALOC) ;
    A = cell2mat(A) ;
    save(['DATAWS/',NAMEDATA,'.mat'],'A');
else
    try
        load(['DATAWS/',NAMEDATA,'.mat'],'A') ;
    end
end



%%%%%%%%%%%%%%%%%%
% Standard SVD
%%%%%%%%%%%%%%%%%%%
NAMECLASS = ['DATAWS/',NAMEDATA,'_class.mat'] ;
if RUN_CLASSIC == 1
    disp('Standard SVD...')
    DATA.RELATIVE_SVD = 1;
    DATA.COMPUTE_V = 0 ;
    tic
    [~,S,~,eSVDstandard]  =SVDT(A,0,DATA) ;
    disp('...End')
    TIME_standard = toc ;
    
    %   TIME_standard = TIME ;
    save(NAMECLASS,'S','TIME_standard','eSVDstandard')
    NORUNC = 0 ;
else
    if exist(NAMECLASS)==0
        NORUNC = 1;
    else
        load(NAMECLASS)
        NORUNC = 0 ;
    end
end

if NORUNC == 0
    
    TRUNCATE = 0;
    
    if  e0>0 & TRUNCATE == 1
        epsilon = e0*sqrt(sum(S.^2)) ;
        SingVsq =  (S.*S) ;
        SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
        normEf2 = sqrt(cumsum(SingVsq)) ; % s_r , s_r +
        T = (sum(normEf2<epsilon)) ;
        K = length(S)-T ;
        eSVDstandard = sqrt(sum(S(K+1:end).^2)) ;
        
        S = S(1:K) ;
        
        
    end
    
    LEGG = 'CLASS.SVD'
    LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVDstandard),' time =',num2str(TIME_standard),' s']
    figure(2)
    hold     on
    h =  plot(log10(S),'r')
    legend(h,LEGG)
    
    
    
    figure(3)
    hold on
    
    SingVsq =  (S.*S) ;
    SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
    %normEf2 = sqrt(sum(SingVsq) - cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
    normEf2 = sqrt(cumsum(SingVsq)) ;
    normEf2 = sort(normEf2,'descend');
    hold on
    lS = log10(normEf2);
    
    h = plot([0:length(S)-1],lS,'r');
    
    xlabel('Modes')
    ylabel('log(error)')
    
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
DATA.COMPUTE_V_SVD = 0 ;
disp('Random. SVD...')
[~,S,~,eSVDrandom]  =RSVDT(A,e0,mu,R,DATA) ;
disp('...End')
TIME = toc
eSVD = eSVDrandom ;
LEGG = ['RSVDT'] ;
LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVD),' time =',num2str(TIME),' s'];
figure(2)
hold     on
h =  plot(log10(S),'b--')
legend(h,LEGG)
xlabel('Modes')
ylabel('log10(S)')

%%%%%%%%%%%%%%%%%%%%%%%%%
legend off
legend show

figure(3)
hold on

SingVsq =  (S.*S) ;
SingVsq = sort(SingVsq);  % s_r, s_{r-1} ... s_1
%normEf2 = sqrt(sum(SingVsq) - cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
normEf2 = sqrt(cumsum(SingVsq)) ;
normEf2 = sort(normEf2,'descend');
hold on
lS = log10(normEf2);

h = plot([0:length(S)-1],lS,'b');

xlabel('Modes')
ylabel('log(error)')

legend(h,LEGG)

if e0 >0
    aaa =  axis;
    h =    plot([aaa(1:2)],log10(e0)*ones(1,2),'g-') ;
    legend(h,['\epsilon'])
end

legend off
legend show
% 
% if RUN_SELECT_COLUMNS ==1
%     [Q B nC c]= RORTH_columns(A,mu,R,DATA) ;
% end