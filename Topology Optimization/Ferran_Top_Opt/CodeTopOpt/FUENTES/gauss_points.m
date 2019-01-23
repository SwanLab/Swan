function [ngaus] = gauss_points(type)

% Cálculo número de puntos de Gauss según tipo elemento

    switch type
        case {'LINEAR_TRIANGLE','LINEAR_TRIANGLE_MIX','LINEAR_TRIANGLE_MIX_COUPLED'}
            ngaus = 1;
        case 'QUAD'
            ngaus = 4;
        otherwise
            error('El elemento no está implementado')
    end
    
end
        