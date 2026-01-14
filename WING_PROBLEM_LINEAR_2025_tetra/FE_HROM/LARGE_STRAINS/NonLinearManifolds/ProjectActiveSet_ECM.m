function [y,dy,d2y] = ProjectActiveSet_ECM(y,cols,VOL,dy,d2y)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/ExtrapWEIGHTS_ECM.pdf
% https://chatgpt.com/share/68c6b960-b340-8013-92d7-a27d7a3d199c
U = y(:,cols) ;
[y(:,cols), PA] = project_simplex_cols(y(:,cols), VOL);
% Optional: map derivatives through the projection’s (piecewise-constant) Jacobian
% (valid as long as the active set is unchanged in a small neighborhood)
dy(:,cols)  = apply_PA(PA, U, dy(:,cols));
d2y(:,cols) = apply_PA(PA, U, d2y(:,cols));