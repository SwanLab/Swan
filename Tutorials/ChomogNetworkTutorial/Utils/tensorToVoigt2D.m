function C_voigt = tensor_to_voigt_2D(C_tensor)
    % Converts a 2x2x2x2 stiffness tensor to a 3x3 Voigt matrix
    %
    % Input:
    %   C_tensor: 4th-order tensor of size (2,2,2,2)
    %
    % Output:
    %   C_voigt: 3x3 Voigt stiffness matrix

    C_voigt = zeros(3,3);
    
    % Normal-normal components
    C_voigt(1,1) = C_tensor(1,1,1,1);   % C11 = C1111
    C_voigt(1,2) = C_tensor(1,1,2,2);   % C12 = C1122
    C_voigt(2,1) = C_tensor(2,2,1,1);   % C21 = C2211
    C_voigt(2,2) = C_tensor(2,2,2,2);   % C22 = C2222
    
    % Normal-shear and shear-normal components
    C_voigt(1,3) = 2 * C_tensor(1,1,1,2);   % C13 = 2*C1112
    C_voigt(2,3) = 2 * C_tensor(2,2,1,2);   % C23 = 2*C2212
    C_voigt(3,1) = 2 * C_tensor(1,2,1,1);   % C31 = 2*C1211
    C_voigt(3,2) = 2 * C_tensor(1,2,2,2);   % C32 = 2*C1222
    
    % Shear-shear component
    C_voigt(3,3) = 4 * C_tensor(1,2,1,2);   % C33 = 4*C1212

end