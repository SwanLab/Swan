function x1 = newton(f,df,x0)
%-------------------------------------------------------------------------%
%      MyBEM 2 - Matthieu Aussal & Francois Alouges - Copyright 2016      %
%                                                                         %
% Ce logiciel MyBEM est la propriete de l'Ecole Polytechnique, tous les   %
% droits lui sont reserves, toute utilisation de ce logiciel est soumise  %
% a l'accord prealable et ecrit de l'Ecole Polytechnique.                 %
%                                                                         %
% This software MyBEM is owned by Ecole polytechnique, all rights are     %
% reserved, any use of this software is subjected to the written, prior   %
% consent of Ecole polytechnique.                                         %
%                                                                         %
% Contact :                                                               %
% matthieu.aussal@polytechnique.edu                                       %
% francois.alouges@polytechnique.edu                                      %
% martin.averseng@polytechnique.edu                                       %
%-------------------------------------------------------------------------%
%
% Creation : 2016.01.01
%
% Last Modification :
% 2016.09.12 
% Synopsis :
% Implementation of Newton method


res  = 1;
while res > 1e-12
    x1 = x0 - f(x0)./df(x0);
    res = max(abs(x1-x0))./max(abs(x0));
    x0 = x1;     
end
x1 = unique(x1);
end