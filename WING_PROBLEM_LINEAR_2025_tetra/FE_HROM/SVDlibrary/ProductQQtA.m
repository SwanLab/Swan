function [QtA   dA  normA] = ProductQQtA(A,Q,alpha,beta,NORESID,TRANSPOSE)

%dbstop('4')
if nargin == 0
    load('tmp1.mat')
    % A = B ;
end


if nargin ==4
    NORESID = 0 ;
    TRANSPOSE = 0 ;
elseif nargin ==5
    TRANSPOSE = 0 ;
end







% [alpha beta] = cellfun(@size,A) ;
% alpha = alpha(:,1) ;
% beta = beta(1,:) ;

if TRANSPOSE == 0
    QtA = 0 ;
else
    QtA = [] ;
end

normA = zeros(size(A)) ;

iacum = 1;


M = size(A,1) ;
N = size(A,2) ;


for i = 1:M
    if TRANSPOSE ==0 ;
        QtAi = [] ;
    else
        QtAi = 0 ;
        iacum = 1;
    end
    for j = 1:N
        
        if ischar(A{i,j})
            SSS = load(B{i,j}) ;
            fff = fieldnames(SSS) ;
            if TRANSPOSE == 0
                alphaLOC_i = size(SSS.(fff{1}),1) ;
                INDICES = iacum:iacum+alphaLOC_i-1 ;
                QtAi = [QtAi   Q(INDICES,:)'*SSS.(fff{1}) ];
            else
                alphaLOC_i = size(SSS.(fff{1}),2) ;
                INDICES = iacum:iacum+alphaLOC_i-1 ;
                QtAi = [QtAi  ; SSS.(fff{1})*Q(:,INDICES) ];
            end
        else
            if TRANSPOSE == 0
                alphaLOC_i = size(A{i,j},1) ;
                INDICES = iacum:iacum+alphaLOC_i-1 ;
                QtAi =  [QtAi   Q(INDICES,:)'*A{i,j}];
            else
                alphaLOC_i = size(A{i,j},2) ;
                INDICES = iacum:iacum+alphaLOC_i-1 ;
                QtAi =  QtAi +  A{i,j}*Q(INDICES,:);
            end
        end
        
        if TRANSPOSE == 1
            iacum = iacum + alphaLOC_i;
        end
        
    end
    
    if TRANSPOSE == 0
        iacum = iacum + alphaLOC_i;
    end
    if TRANSPOSE == 0
        QtA = QtA + QtAi ;
    else
        QtA = [QtA ; QtAi] ;
    end
end



dA = cell(size(A)) ;
iacumROW = 1;
normA = [] ;
if NORESID == 0
    for i = 1:size(A,1)
        iacumCOL = 1 ;
        INDICESrow = iacumROW:iacumROW+alpha(i)-1 ;
        for j = 1:size(A,2)
            INDICEScol = iacumCOL:iacumCOL+beta(j)-1 ;
            if ischar(A{i,j})
                SSS = load(B{i,j}) ;
                fff = fieldnames(SSS) ;
                %      if NORESID ==0
                dAij = SSS.(fff{1}) - Q(INDICESrow,:)*QtA(:,INDICEScol) ;
                %     else
                %        dAij =   QtA(:,INDICEScol)'*Q(INDICESrow,:)' ;
                %   end
                
                normA(i,j) = norm(SSS.(fff{1}),'fro') ;
                error('Completar esta parte')
                
                
                
            else
                
                %   if NORESID ==0
                dA{i,j} = A{i,j} -  Q(INDICESrow,:)*QtA(:,INDICEScol) ;
                %  else
                %     dA{i,j} = QtA(:,INDICEScol)'*Q(INDICESrow,:)' ;
                % end
                normA(i,j) = norm(A{i,j},'fro') ;
                
            end
            iacumCOL = iacumCOL + beta(j);
        end
        iacumROW = iacumROW + alpha(i);
    end
    
    % if NORESID == 1
    %
    %     dA = dA' ;
    % end
    
    
end
