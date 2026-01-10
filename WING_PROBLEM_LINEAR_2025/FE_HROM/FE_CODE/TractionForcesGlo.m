function Ftrac = TractionForcesGlo(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb,DATA)

if nargin == 0
    load('tmp.mat')
end

if DATA.VECTcode  == 0
    % Standard way (elementwise)
    Ftrac = FtracCOMP(COOR,CNb,TypeElementB,Fpnt,Tnod);
else
    % Vectorized code
    Ftrac = FtracCOMPvect(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb);
end

if exist(DATA.nameWORKSPACE,'file')==0
    APPEND = '' ;
else
    APPEND = '-append' ;
end

if DATA.STORE_STIFFNESS ==1
    save(DATA.nameWORKSPACE,'Ftrac','CNb','Fpnt','Tnod',APPEND);
elseif DATA.STORE_STIFFNESS ==2   %&& DATA.DO_NOT_STORE_GEOMETRIC_INFORMATION == 0
    save(DATA.nameWORKSPACE,'Ftrac','CNb','TypeElementB','Fpnt','Tnod','CONNECTb',APPEND)
end