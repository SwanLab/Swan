function      PERMUT = PermutationConnHexahedra27
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/01_8nodeHEX.mlx
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/03_20nodeHEX.mlx
% Adaptation of /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/PermutationConnHexahedra8.m
% BASIC PERMUTATIONS
% -------------------
PERMUT_base = cell(1,8) ; 
% PERMUTATIONS CONSIDERING THAT THE BOTTOM REFERENCE PLANE IN THE X3MIN
% PLANE
PERMUT_base{1} =  [1,2,3,4,   5,6,7,8,  9,10,11,12,  13,14,15,16,  17,18,19,20,   21,  22,23,24,25,  26,   27] ;
PERMUT_base{2} =  [2,3,4,1    6,7,8,5   10,11,12,9   14,15,16,13   18,19,20,17,   21,  23,24,25,22,  26    27] ;                       
PERMUT_base{3} =  [3,4,1,2    7,8,5,6   11,12,9,10   15,16,13,14   19,20,17,18,   21,  24,25,22,23,  26   27 ] ; 
PERMUT_base{4} =  [4,1,2,3    8,5,6,7   12,9,10,11   16,13,14,15   20,17,18,19,   21,  25,22,23,24,  26,  27  ] ; 
% PERMUATIONS CONSIDERING THAT THE BOTTOM REFERENCE PLANE IS THE X3MAX
% PLANE
PERMUT_base{5} = [8,7,6,5,   4,3,2,1,  19,18,17,20,  16,15,14,13,   11,10,9,12,   26,  24,23,22,25,   21, 27] ;  
PERMUT_base{6} = [7,6,5,8    3,2,1,4   18,17,20,19   15,14,13,16    10,9,12,11,    26,  23,22,25,24,   21, 27] ;
PERMUT_base{7} = [6,5,8,7    2,1,4,3   17,20,19,18   14,13,16,15    9,12,11,10,   26,  22,25,24,23,   21, 27] ; 
PERMUT_base{8} = [5,8,7,6    1,4,3,2   20,19,18,17   13,16,15,14    12,11,10,9,    26,  25,24,23,22   21, 27] ;  

% INITIAL CONFIGURATIONS (depending on the plane one chooses as "bottom" plane, x1min, x2min or x3min)
% -----------------------------
INIT_CONFIG = cell(1,3) ; 
INIT_CONFIG{1} = 1:27 ;    % X3MIN 
INIT_CONFIG{2} = [2,1,5,6,  3,4,8,7,   9,13,17,14,   10,12,20,18,  11,16,19,15,  ...
                  22,21,25,26,23,24,27 ] ; % X2MIN

INIT_CONFIG{3} = [2,6,7,3,  1,5,8,4,  14,18,15,10   9,17,19,11,   13,20,16,12,...
                  23,22,26,24,21,25,27 ] ; % X1MAX

% P21 - 1234 -->2673 --> 23
% P22 - 1256 -->2615 --> 22
% P23 -- 2367 ---6758 -- 26

PERMUT= cell(1,24) ; 

iacum = 1 ; 

for iconf = 1:length(INIT_CONFIG)
    CCONF = INIT_CONFIG{iconf} ; 
    for iperm = 1:length(PERMUT_base)
        PERMUT{iacum} = CCONF(PERMUT_base{iperm}) ; 
        iacum = iacum + 1; 
    end
end


 
