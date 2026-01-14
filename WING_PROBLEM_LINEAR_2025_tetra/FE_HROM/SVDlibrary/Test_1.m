clc
clear all

%%
GENERATE_AGAIN = 1 ;
R= 0 ;
e0 = 0.2;
muREL = 0.2 ;
DATA.RELATIVE_SVD = 1;
RANK_PRE = 300 ;
SIZE_A =0.3;     % Size matrix (Gbbytes)
p =1;   % Number of partitions along rows
q =  5 ;  % Number of partitions along columns
ratioNM  =10 ;
N = ceil(sqrt(SIZE_A/8*1e9)/sqrt(ratioNM))  ; % Number of columns
M = ceil( SIZE_A/8*1e9/N) ;  % Number of rows
DATASTORE = ['DATAWS/DATA_Amatrix/'] ;
NAMEWS_ALL = ['DATAWS/A_ws.mat']

beta = ceil(N/q) ;
beta = repmat(beta,1,q) ;
alpha = M ;
MAXSTORE = 4.1;

SEM = randn(alpha(1),RANK_PRE) ;

METHODS = {'CLASSICAL_SVD','RSVD'}
COLORS = {'r','b--'}

figure(1)
hold     on
xlabel('Modes')
ylabel('log10(S)')

for  imethod = 1:length(METHODS)
    if imethod > 1
        
        GENERATE_AGAIN = 0 ;
    end
    
    
    if   SIZE_A > MAXSTORE  ||  (SIZE_A < MAXSTORE  & GENERATE_AGAIN==1)
        disp('Generation')
        A = cell(p,q) ;
        for i = 1:length(alpha)
            disp(['i=',num2str(i)])
            for j = 1:length(beta) ;
                if  SIZE_A > MAXSTORE
                    disp(['j=',num2str(j)])
                    nameWS =[DATASTORE,'A_',num2str(i),'_',num2str(j),'.mat'] ;
                    
                    if GENERATE_AGAIN == 1
                        Aij = SEM*rand(size(SEM,2),beta(i)) ;
                        disp('Saving...')
                        save(nameWS,'Aij') ;
                        disp('End...')
                    end
                    
                    A{i,j} = nameWS ;
                    
                    
                else
                    disp(['j=',num2str(j)])
                    A{i,j} = SEM*rand(size(SEM,2),beta(i)) ;
                end
            end
        end
        
        disp('End Generation')
        
        disp('SAving ...')
        save(NAMEWS_ALL,'A')
        disp('End')
        
    else
        disp('Loading ...')
        load(NAMEWS_ALL,'A')
        disp('End ...')
    end
    
    
    switch METHODS{imethod }
        case 'CLASSICAL_SVD'
            disp('Classic SVD')
            A =cell2mat(A) ;
            tic
            [U,S,V,eSVDclassic]  =SVDT(A,e0,DATA) ;
            TIME = toc ;
            LEGG = 'CLASS.SVD'
            eSVD = eSVDclassic ;
        case 'RSVD'
            A = cell2mat(A);
            a = norm(A,'fro') ;
            if muREL >0
                mu = muREL*a ;
            else
                mu = 0 ;
            end
            tic
            
            
            e0 = e0*norm(A,'fro') ;
            DATA.RELATIVE_SVD =0 ;
            [U,S,V,eSVDrandom]  =RSVDT(A,e0,mu,R,DATA) ;
            TIME = toc
            eSVD = eSVDrandom ;
            LEGG = ['RSVDT'] ;
    end
    LEGG = [LEGG,' rank=',num2str(length(S)),' error = ',num2str(eSVD),' time =',num2str(TIME),' s']
    figure(1)
    hold     on
    
    h =  plot(log10(S),COLORS{imethod})
    legend(h,LEGG)
    
    
end

legend off
legend show
