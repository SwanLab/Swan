function [U,S,V,DATAOUT] = SVDcol(B,NMODES)
% DATA: Partitioned matrix B = [B1 B2 ... Bq],  }),
% and NMODES 
% RESULT: Orthogonal matrices U,V, and S = diag(s1,s2,,,sk)
% 5-Oct-2016, Joaquín A. Hernández, UPC-CIMNE
% Adapted 21-Nov-2017
% -------------------------------------------------------
 if nargin == 0
     aaa = rand(8) ; 
    B = [aaa,5*aaa] ; 
    NMODES  = [2 2]  ; 
 end
 TIME = tic ; 
 % ---------------------------------------------
 [L,G] = SVDloop(B,beta,epsilon) ; % SVD of all submatrices
 % ----------------------------------------------
 TIME = toc(TIME) ; 
 disp(['TIME = ',num2str(TIME)]) ; DATAOUT.TIMEloop = TIME ; 
 disp('----------------------------------')
 disp('Second step')
 disp(['Size L =',num2str(size(L,1)),' x ',num2str(size(L,2)),'  (of ',num2str(N),' rows)'])
 TIME = tic ; 
 % ----------------------------------
 % --------------------------------------
 [U,S,Vbar] = SVD(L,0) ; 
 V = G*Vbar ; 
 % -----------------------------------------
 % -----------------------------------------------
 TIME = toc(TIME) ; 
 disp(['TIME = ',num2str(TIME)]) ; DATAOUT.TIME2step = TIME ; 
 disp('----------------------------------')