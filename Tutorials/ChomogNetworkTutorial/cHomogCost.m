function [J, grad] = cHomogCost(opt_cHomog, type, x)

    % Fetch output and jacobian of the network
    Y = opt_cHomog.computeOutputValues(x');
    dY = opt_cHomog.computeGradient(x');
    
    % Fetch the Chomog component to optimize
    C_11 = Y(1);
    C_12 = Y(2);
    C_22 = Y(3);
    C_33 = Y(4);
    dC_11 = dY(:, 1)';
    dC_12 = dY(:, 2)';
    dC_22 = dY(:, 3)';
    dC_33 = dY(:, 4)';

    % Calculate the cost and gradient
    switch type
        case 'maxHorzStiffness'
            J = 1 / C_11;
            grad = - 1 / C_11^2 .* dC_11;
        case 'maxVertStiffness'
            J = 1 / C_22;
            grad = - 1 / C_22^2 .* dC_22;
        case 'maxShearStiffness'
            J = 1 / C_33;
            grad = - 1 / C_33^2 .* dC_33;
        case 'maxIsoStiffness'
            J = 1 / (C_11 + C_22);
            grad = - (dC_11 + dC_22) / (C_11 + C_22)^2;
        case 'minHorzStiffness'
            J = C_11;
            grad = dC_11;
        case 'minVertStiffness'
            J = C_22;
            grad = dC_22;
        case 'minShearStiffness'
            J = C_33;
            grad = dC_33;
        case 'minIsoStiffness'
            J = C_11 + C_22;
            grad = dC_11 + dC_22;
        case 'maxAuxetic'
            J = C_12 / C_11;
            grad = dC_12 / C_11 - C_12 / C_11^2 * dC_11;
    end

end