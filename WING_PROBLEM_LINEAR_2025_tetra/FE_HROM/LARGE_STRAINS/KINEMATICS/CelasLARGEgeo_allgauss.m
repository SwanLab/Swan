function celasLARGEgeo = CelasLARGEgeo_allgauss(StwoST,ndim)
% Assembly Geometric celasLARGE Matrix
% See % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 18
% See also /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/Kgeometric.m
%
if nargin == 0
    load('tmp.mat')
end

if ndim == 2
    nF = 4 ;
    nstrain = 3;
    nelem_ngaus = length(StwoST)/nstrain  ;
    celasLARGEgeo = zeros(nF*nelem_ngaus,nF)  ;
    
    CROWS = cell(1,nF) ;
    for icols  =1:nF
        CROWS{icols} = icols:nF:size(celasLARGEgeo,1) ;
    end
    
    SROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        SROWS{icols} = icols:nstrain:size(StwoST,1) ;
    end
    
    %     % celasGEO =
    %
    % [S1,  0, S3,  0]
    % [ 0, S2,  0, S3]
    % [S3,  0, S2,  0]
    % [ 0, S3,  0, S1]
    
    % ROW 1 =  [S1,  0, S3,  0]
    celasLARGEgeo(CROWS{1},1) = StwoST(SROWS{1}) ;
    celasLARGEgeo(CROWS{1},3) = StwoST(SROWS{3}) ;
    % ROW 2 = [ 0, S2,  0, S3]
    celasLARGEgeo(CROWS{2},2) = StwoST(SROWS{2}) ;
    celasLARGEgeo(CROWS{2},4) = StwoST(SROWS{3}) ;
    % ROW 3 =[S3,  0, S2,  0]
    celasLARGEgeo(CROWS{3},1) = StwoST(SROWS{3}) ;
    celasLARGEgeo(CROWS{3},3) = StwoST(SROWS{2}) ;
    % ROW 4 =[ 0, S3,  0, S1]
    celasLARGEgeo(CROWS{4},2) = StwoST(SROWS{3}) ;
    celasLARGEgeo(CROWS{4},4) = StwoST(SROWS{1}) ;
    
    
else
    nF = 9 ;
    nstrain = 6;
    nelem_ngaus = length(StwoST)/nstrain  ;
    celasLARGEgeo = zeros(nF*nelem_ngaus,nF)  ;
    
    CROWS = cell(1,nF) ;
    for icols  =1:nF
        CROWS{icols} = icols:nF:size(celasLARGEgeo,1) ;
    end
    
    SROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        SROWS{icols} = icols:nstrain:size(StwoST,1) ;
    end
    
    %
    %     celasGEO =
    %
    % [S1,  0,  0,  0, S5, S6,  0,  0,  0]  1
    % [ 0, S2,  0, S4,  0,  0,  0,  0, S6]  2
    % [ 0,  0, S3,  0,  0,  0, S4, S5,  0]  3
    % [ 0, S4,  0, S3,  0,  0,  0,  0, S5]  4
    % [S5,  0,  0,  0, S3, S4,  0,  0,  0]  5
    % [S6,  0,  0,  0, S4, S2,  0,  0,  0]  6
    % [ 0,  0, S4,  0,  0,  0, S2, S6,  0]  7
    % [ 0,  0, S5,  0,  0,  0, S6, S1,  0]  8
    % [ 0, S6,  0, S5,  0,  0,  0,  0, S1]  9
    
    
    % ROW 1 =  [S1,  0,  0,  0, S5, S6,  0,  0,  0]
    celasLARGEgeo(CROWS{1},1) = StwoST(SROWS{1}) ;
    celasLARGEgeo(CROWS{1},5) = StwoST(SROWS{5}) ;
    celasLARGEgeo(CROWS{1},6) = StwoST(SROWS{6}) ;
    % ROW 2 = [ 0, S2,  0, S4,  0,  0,  0,  0, S6]
    celasLARGEgeo(CROWS{2},2) = StwoST(SROWS{2}) ;
    celasLARGEgeo(CROWS{2},4) = StwoST(SROWS{4}) ;
    celasLARGEgeo(CROWS{2},9) = StwoST(SROWS{6}) ;
    % ROW 3 =[ 0,  0, S3,  0,  0,  0, S4, S5,  0]
    celasLARGEgeo(CROWS{3},3) = StwoST(SROWS{3}) ;
    celasLARGEgeo(CROWS{3},7) = StwoST(SROWS{4}) ;
    celasLARGEgeo(CROWS{3},8) = StwoST(SROWS{5}) ;
    % ROW 4 =[ 0, S4,  0, S3,  0,  0,  0,  0, S5] 
    celasLARGEgeo(CROWS{4},2) = StwoST(SROWS{4}) ;
    celasLARGEgeo(CROWS{4},4) = StwoST(SROWS{3}) ;
    celasLARGEgeo(CROWS{4},9) = StwoST(SROWS{5}) ;
    % ROW 5 = [S5,  0,  0,  0, S3, S4,  0,  0,  0]
    celasLARGEgeo(CROWS{5},1) = StwoST(SROWS{5}) ;
    celasLARGEgeo(CROWS{5},5) = StwoST(SROWS{3}) ;
    celasLARGEgeo(CROWS{5},6) = StwoST(SROWS{4}) ;
    % ROW 6 = [S6,  0,  0,  0, S4, S2,  0,  0,  0]
    celasLARGEgeo(CROWS{6},1) = StwoST(SROWS{6}) ;
    celasLARGEgeo(CROWS{6},5) = StwoST(SROWS{4}) ;
    celasLARGEgeo(CROWS{6},6) = StwoST(SROWS{2}) ;
    % ROW 7 = [ 0,  0, S4,  0,  0,  0, S2, S6,  0]
    celasLARGEgeo(CROWS{7},3) = StwoST(SROWS{4}) ;
    celasLARGEgeo(CROWS{7},7) = StwoST(SROWS{2}) ;
    celasLARGEgeo(CROWS{7},8) = StwoST(SROWS{6}) ;
    % ROW 8 = [ 0,  0, S5,  0,  0,  0, S6, S1,  0]
    celasLARGEgeo(CROWS{8},3) = StwoST(SROWS{5}) ;
    celasLARGEgeo(CROWS{8},7) = StwoST(SROWS{6}) ;
    celasLARGEgeo(CROWS{8},8) = StwoST(SROWS{1}) ;
    % ROW 9 = [ 0, S6,  0, S5,  0,  0,  0,  0, S1]
    celasLARGEgeo(CROWS{9},2) = StwoST(SROWS{6}) ;
    celasLARGEgeo(CROWS{9},4) = StwoST(SROWS{5}) ;
    celasLARGEgeo(CROWS{9},9) = StwoST(SROWS{1}) ;
    
    
    
end