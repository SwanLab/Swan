function [U,S,V,ETIME,RankR,e_svd] = RSVDqpCOMP(A,epsilon,DATA)

if nargin == 0
    load('tmp1.mat')
end
% First call to RSVDqp
ETIME.Al = [] ;
ETIME.Ar = [] ;
U = [] ; S = [] ; V = [] ;

RankR = [] ;


[alpha beta] = DetermineAlphaBeta(A) ;


EPSILON_GLO = DATA.EPSILON_GLO ;
DATA.EPSILON_GLO = 0 ;
disp('..........................')
disp('Constrained SVD')
disp('..........................')

DATA.alphaA = alpha ;
DATA.betaA = beta ;

DATA = DefaultField(DATA,'DoNotSepateROWCOLS',0) ;
DATA = DefaultField(DATA,'COLR',[]) ;
DATA = DefaultField(DATA,'ROWR',[]) ;
DATA = DefaultField(DATA,'ProjectONTO_rangeA',0) ;


DATA = DefaultField(DATA,'ProjectQr_onto_colA',0) ;

if  DATA.DoNotSepateROWCOLS == 1
    
    %    error('This method IS BY NO MEANS EFFICIENT. ERASE IT !!!')
    DATA.NoCalculateQTAWhenEPS1 = 0 ;
    [Ur,Sr,Vr] = RSVDqp(A,epsilon,DATA)  ;
    
    Qr = Ur ;
    
    
    if DATA.ProjectQr_onto_colA == 1
        epsilonLOC = zeros(size(epsilon)) ;
        [U,~,~,~,~,~,DATAOUT]  = RSVDqp(A,epsilonLOC,DATA)  ;
        
        Qr = U*(U'*Qr) ;
        Qr = ORTH(Qr,DATAOUT.mu) ;
        
    end
    
    
    
else
    % Constrained columns
    % -----------------------------------------
    
    
    
    DATA.SetGammaToZero = 1  ;
    if ~isempty(DATA.COLR)
        epsilonLOC = ones(size(epsilon));
        epsilonLOC(:,DATA.COLR) =epsilon(:,DATA.COLR);
    else
        SSS = sum(epsilon,1) ;
        COLexact =  find(SSS ==0);
        epsilonLOC = ones(size(epsilon)) ;
        epsilonLOC(:,COLexact) = 0;
    end
    % ---------------------------------------------------------
    Qr = [] ;
    if ~all(all(epsilonLOC==1))
        DATA.NoCalculateQTAWhenEPS1 = 1 ;
        [Qr] = RSVDqp(A,epsilonLOC,DATA)  ;
    end
    % Constrained rows
    % -----------------------------------------
    if ~isempty(DATA.COLR)
        epsilonLOC = ones(size(epsilon));
        epsilonLOC(DATA.ROWR,:) =epsilon(DATA.ROWR,:);
    else
        SSS = sum(epsilon,2) ;
        ROWexact =  find(SSS ==0);
        epsilonLOC = ones(size(epsilon)) ;
        epsilonLOC(ROWexact,:) = 0;
    end
    Qrr = [] ;
    if ~all(all(epsilonLOC==1))
        DATA.NoCalculateQTAWhenEPS1 =1 ;
        [Qrr,~,~,~,~,~,DATAOUT] = RSVDqp(A,epsilonLOC,DATA)  ;
    end
    %%%
    
    
    
    
    
    if DATA.ProjectQr_onto_colA == 1
        if ~isempty(Qrr)
            if ~isempty(Qr)
                %  Qrr = ORTH(Qrr- Qr*(Qr'*Qrr),DATAOUT.mu) ;
                Qr = [Qr Qrr] ;
            else
                Qr = Qrr ;
            end
        end
        
        epsilonLOC = zeros(size(epsilon)) ;
        [U,~,~,~,~,~,DATAOUT]  = RSVDqp(A,epsilonLOC,DATA)  ;
        
        Qr = U*(U'*Qr) ;
        Qr = ORTH(Qr,DATAOUT.mu) ;
        
    else
        if ~isempty(Qrr)
            if ~isempty(Qr)
                Qrr = ORTH(Qrr- Qr*(Qr'*Qrr),DATAOUT.mu) ;
                Qr = [Qr Qrr] ;
            else
                Qr = Qrr ;
            end
        end
        
        
    end
    
end


NORESID  = 1 ;
TRANSPONSE = 0;
[QrTA] = ProductQQtA(A,Qr,alpha,beta,NORESID,TRANSPONSE) ;
[Ur Sr Vr] = RSVDT(QrTA,0) ;
Ur = Qr*Ur ;


%
% end



NORESID  = 1 ;
TRANSPONSE = 0;
[QrTA,dA,normAoriginal] = ProductQQtA(A,Ur,alpha,beta);

nA = norm(normAoriginal,'fro') ; % Norm of the initial residual
mu =   (max([sum(alpha),sum(beta)])*eps(nA))  ;  ; % MAchine precision parameter

epsilonLOC = zeros(size(epsilon)) ;


DATA.normAoriginal = normAoriginal ;

[Ul,Sl,Vl,ETIME.Al] = RSVDqp(dA,epsilonLOC,DATA)  ;
% else
%
%     AA = cell2mat(A) ;
%
%     dA = AA - Ur*(Ur'*AA) ;
%     normAoriginal  =norm(AA,'fro') ;
%     [Ul,Sl,Vl ] = RSVDqp(dA,0)  ;
%     [Ur,Sr,Vr ] = RSVDqp(Ur*(Ur'*AA),0)  ;
% end

% Computing Al




U = [Ur Ul] ;
S = [Sr; Sl] ;
V = [Vr Vl] ;






DATA.EPSILON_GLO = EPSILON_GLO ;

c = norm(normAoriginal,'fro') ;
%
% NSQ = cumsum(S.^2) ;
% LLL = find(NSQ >c^2) ;
% R =length(S)-length(LLL) ;
% U = U(:,1:R);
% V = V(:,1:R) ;
% S = S(1:R) ;

%dbstop('159')
mu = (max([sum(alpha) sum(beta)])*eps(c))  ;  % Machine epsilon parameter
eLOC = DATA.EPSILON_GLO ;
if  eLOC == 0
    e0 = mu ;
else
    e0 = c*eLOC ;
end

% -----------------------------------------
% Estimation rank A

if  e0<=mu
    R = length(find(S>mu)) ;
else
    %  disp(['Rank matrix =',num2str(RankMatrix)]) ;
    SingVsq =  (S.*S) ;
    SingVsq = SingVsq(length(SingVsq):-1:1);  % s_r, s_{r-1} ... s_1
    normEf2 = sqrt(cumsum(SingVsq)) ; % s_r , s_r + s_{r-1} ... s_r +s_{r-1}+ s_{r-2} ... + s_1
    tol = e0 ;
    if tol<=mu
        R = length(S) ;
    else
        T = (sum(normEf2<tol)) ;
        R = length(S)-T ;
    end
end
% Actual error
e_svd = sqrt(sum(S(R+1:end).^2)) ;
U = U(:,1:R);
S = S(1:R) ;
if ~isempty(V)
    V = V(:,1:R) ;
end
disp(['Truncated rank=',num2str(R)]) ;



SingVsq =  (S.*S) ;
DENOM = sum(SingVsq) ;    NOMINAD = DENOM-cumsum(SingVsq) ;
ERRORsvd = sqrt(NOMINAD/DENOM);


 DATA =DefaultField(DATA,'PLOT_GRAPH_SVD',0) ; 

 if DATA.PLOT_GRAPH_SVD == 1
figure(100)
hold on
xlabel('Mode')
ylabel('Error SVD (%)')

plot(ERRORsvd*100,'b')

title(['Rank Restr Matrix = ',num2str(length(Sr))])

aaa = axis ;

r =length(Sr) ;

plot([r r],[aaa(3:4)],'r-')
 end

RankR = length(Sr)  ;

% U = Ur;
% V = Vr;
% S = Sr;

