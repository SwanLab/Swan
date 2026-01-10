function [L,G,normBapprox] = RSVDloop(B,NMODES,DATA)
% DATA: Partitioned matrix B = {B1 B2 ... Bq},
 % RESULT: L = [L1 L2 ... Lq], and G = diag(G1,G2 ...Gq), where Li = Di*Hi
% and Bi = Di*Hi*Di^T  
if nargin == 0
    load('tmp1.mat')
elseif nargin == 2
    DATA = [] ; 
end

if iscell(B)
DATA= DefaultField(DATA,'TOL_BLOCK_SVD',cell(size(B))) ; 
end

disp('----------------')
disp('LOOP OVER BLOCKS')
disp('----------------')
q = length(B) ;
L = cell(1,q) ;
G = cell(1,q) ;


if ~iscell(B)
    M = size(B,1) ;
    B = mat2cell(B,M,beta) ;
end

normBapprox = zeros(1,q) ;

for i=1:q
    disp('------------------------------------------')
    disp(['i = ',num2str(i), '  of ',num2str(q)])
    TTT =tic ;
    if   isempty(NMODES{i})
        estimationRANK = 0 ;
    else
        Rsup = min(size(B{i})) ;
        estimationRANK = min(Rsup,NMODES{i}) ;
    end
    if isnumeric(B{i})
        TOL_LOC = DATA.TOL_BLOCK_SVD{i} ; 
        [Di,Hi,Gi] = RSVDT(B{i},TOL_LOC,[],estimationRANK,DATA); 
    elseif iscell(B{i})
        
        % Invoking partitioned method, with tolerance given by DATA.TOL_BLOCK_SVD
 
        TOL_BLOCK = DATA.TOL_BLOCK_SVD{i} ; 
         [Di,Hi,Gi,ETIME,eSVD,RankMatrix ] = RSVDqp(B{i},TOL_BLOCK)  ; 
       
    else
        disp('Retrieving from memory ...')
        load(B{i}) ; % Retrieved from memory
        disp('Done')
         [Di,Hi,Gi] = RSVDT(B{i},0,[],estimationRANK); 
     
        disp('Deleting...')
        delete(B{i});
        
    end
    
    if ~isempty(NMODES{i})
    NMODES_i = min(NMODES{i},length(Hi));
    else
        NMODES_i  = length(Hi); 
    end
    
    Di = Di(:,1:NMODES_i );
    Hi = Hi(1:NMODES_i) ;
    Gi = Gi(:,1:NMODES_i) ;
    normBapprox(i) = sqrt(sum(Hi.^2)) ;
    %    dbstop('40')
    Ni = size(B{i},2) ;
    
    TTT= toc(TTT);
    
    L{i} = bsxfun(@times,Di',Hi)' ;
    %  dbstop('13')
    G{i} = sparse(Gi) ;
    disp(['K = ',num2str(length(Hi)),' of ',num2str(Ni),' columns'])
    disp(['Time = ',num2str(TTT)])
end

% Writing G in sparse format
% G = diagonal(G{1},G{2},...G{q})
%dbstop('26')
L = cell2mat(L) ;
G = blkdiag(G{:}) ;
