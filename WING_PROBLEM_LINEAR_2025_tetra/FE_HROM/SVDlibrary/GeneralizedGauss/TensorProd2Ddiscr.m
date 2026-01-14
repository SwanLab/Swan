function [MPOINTS Mglo x w]  =TensorProd2Ddiscr(M,xLIM)

if nargin == 0
    M = [3,3]; 
    xLIM = [-1 1 ; -1 1] ; 
elseif    nargin ==1
    xLIM = [-1 1 ; -1 1]; 
    
end

%     %
%     if M(1) ~=M(2)
% %         error('Not valid discretization  for tensor product rules ')
%     end
%     doffXLIM = xLIM-[-1 1 ; -1 1];
%     if norm(doffXLIM)~=0
%         error('Not valid limits (it should be a biunit square)')
%     end
    
    
    [x1d, w1d] = GaussQuad(M(1), xLIM(1,1), xLIM(1,2)) ;
    [x2d, w2d] = GaussQuad(M(2), xLIM(2,1), xLIM(2,2)) ;
    %
    %xxx =[]  ;
    %yyy = [] ;
    iacum = 1;
    w = [] ;
    for i = 1:length(w1d)
        for j = 1:length(w2d)
            %       xxx(iacum)  = x1d(i) ;
            %      yyy(iacum) =  x1d(j) ;
            w(iacum) = w1d(i)*w2d(j) ;
            iacum = iacum+1 ;
        end
    end
    
    x{1} = x1d' ;
    x{2} = x2d' ;
    MPOINTS = M ;
    Mglo = MPOINTS(1)*MPOINTS(2) ;
    w = w' ;
    