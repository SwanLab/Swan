function DATAFITTING = GetCoefficientesFitting(xINI,G,DATA)  

ndim = size(xINI,2) ; 

if ndim ==1 
     DATAFITTING = GetCoefficientesFitting1D(xINI,G) ; 
elseif ndim ==2 
     DATAFITTING = GetCoefficientesFitting2D(xINI,G,DATA) ; 
elseif ndim ==3
    DATAFITTING = GetCoefficientesFitting3D(xINI,G,DATA) ; 
else
    error('Option not implemented')
end