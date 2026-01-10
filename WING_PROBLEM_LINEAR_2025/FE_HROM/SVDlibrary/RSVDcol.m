function [U,S,V,L] = RSVDcol(B,NMODES,DATA)
% DATA: Partitioned matrix B = [B1 B2 ... Bq],  }),
% and NMODES 
% RESULT: Orthogonal matrices U,V, and S = diag(s1,s2,,,sk)
% 5-Oct-2016, Joaquín A. Hernández, UPC-CIMNE
% Adapted 21-Nov-2017
% -------------------------------------------------------
 if nargin == 0
     load('tmp1.mat')
 elseif nargin == 2 
     DATA = [] ; 
 end
  % ---------------------------------------------
 [L,G] = RSVDloop(B,NMODES,DATA) ; % SVD of all submatrices
 % ----------------------------------------------
 %  disp(['TIME = ',num2str(TIME)]) ; DATAOUT.TIMEloop = TIME ; 
%  disp('----------------------------------')
%  disp('Second step')
% disp(['Size L =',num2str(size(L,1)),' x ',num2str(size(L,2)),'  (of ',num2str(N),' rows)'])
 
 % ----------------------------------
 % --------------------------------------
 DATA = DefaultField(DATA,'TOL_SVD_GLO',0) ;
 DATALOC.RELATIVE_SVD = 1;
 
 [U,S,Vbar] = RSVDT(L,DATA.TOL_SVD_GLO,0,0,DATALOC) ;
 V = G*Vbar ;
 % -----------------------------------------
 % -----------------------------------------------
 
%  disp(['TIME = ',num2str(TIME)]) ; DATAOUT.TIME2step = TIME ; 
%  disp('----------------------------------')