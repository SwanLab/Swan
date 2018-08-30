function gradient = ritz_representation_of_gradient(gradient,epsilon_P0_gradient,Msmooth,Stiff_smooth,algorithm,ritz_gradient_representation_case)

switch algorithm
    case 'gradient'
        switch ritz_gradient_representation_case
            case 'L2'
            gradient = gradient;
        case 'H1'
            gradient = (epsilon_P0_gradient^2*Stiff_smooth + Msmooth)\(Msmooth*gradient);
        end
end

end