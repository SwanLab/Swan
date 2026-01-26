function [dimREF signREF] = Faces2Index(FACE)

switch  FACE
    case 'xMIN'
        dimREF = 1;
        signREF = 1 ;
    case 'xMAX'
        dimREF = 1;
        signREF = -1 ;
    case 'yMIN'
        dimREF = 2;
        signREF = 1 ;
    case 'yMAX'
        dimREF = 2;
        signREF = -1 ;
    case 'zMIN'
        dimREF = 3;
        signREF = 1 ;
    case 'zMAX'
        dimREF = 3;
        signREF = -1 ;
end