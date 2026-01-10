function [L,G,normBapprox] = SVDloop(B,beta,epsilon)
% DATA: Partitioned matrix B = {B1 B2 ... Bq},
% error threshold epsilon (1 x q)
% RESULT: L = [L1 L2 ... Lq], and G = diag(G1,G2 ...Gq), where Li = Di*Hi
% and Bi = Di*Hi*Di^T + ERROR(epsilon)
if nargin == 0
    load('tmp.mat')
end
disp('----------------')
disp('LOOP OVER BLOCKS')
disp('----------------')
q = length(beta) ;
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
    if isnumeric(B{i})
        
        
        %   dbstop('32')
        [Di,Hi,Gi,DATAOUT] = SVD(B{i},epsilon(i)); 
        normBapprox(i) = sqrt(sum(Hi.^2)) ; 
        %    dbstop('40')
        Ni = size(B{i},2) ;
    else
        disp('Retrieving from memory ...')
        load(B{i}) ; % Retrieved from memory
        disp('Done')
        [Di,Hi,Gi] = SVD(Bi,epsilon(i));
        Ni = size(Bi,2) ;
        % Delete from memory
        disp('Deleting...')
        delete(B{i});
        
    end
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
