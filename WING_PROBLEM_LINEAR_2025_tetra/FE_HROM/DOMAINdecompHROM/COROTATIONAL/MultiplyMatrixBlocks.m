function C = MultiplyMatrixBlocks(A,B)
% Multiplication in a blockwise fashion 
% C = A*B, where 
% A = [A1;A2; ...], B = [B1;B2; ...]
% Ci = Ai*Bi ; 
% JAHO, 3-Dec-2024, UPC, Terrassa 
% ---------------------------------
ndim = size(A,2) ; 
  C = zeros(size(B)) ; 
     nrows = size(C,1) ; 
     for iiiLOC = 1:ndim
         iii = iiiLOC:ndim:nrows ;
         for jjj = 1:ndim
             for kkkLOC = 1:ndim
                 kkk = kkkLOC:ndim:nrows;
                 C(iii,jjj) =   C(iii,jjj)  + A(iii,kkkLOC).*B(kkk,jjj) ;
             end
         end
     end