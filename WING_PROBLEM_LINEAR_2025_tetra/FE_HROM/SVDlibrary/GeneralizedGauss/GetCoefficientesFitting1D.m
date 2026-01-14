function DATAFITTING = GetCoefficientesFitting1D(xINI,G)  
if nargin == 0
    load('tmp.mat')
end
nmodes= size(G,1) ; % Number of functions to be integrated 
DATAFITTING.spline_G = cell(1,nmodes) ; % Coefficients splines for each column
DATAFITTING.spline_derG = cell(1,nmodes) ;          % Coefficients derivatives splines for each column
 for i = 1:nmodes
    DATAFITTING.spline_G{i} = spline(xINI,G(i,:)) ;
    DATAFITTING.spline_derG{i} =  fnder( DATAFITTING.spline_G{i}  ) ;
end