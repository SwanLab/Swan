function [fobj,gradient] = assamble_cost_gradient(cost_mechanical,gradient_mechanical,constraint_Perimeter,gradient_constraint_Perimeter,...
                                                    lambda_Perimeter,penalty_Perimeter)

cost_Perimeter = lambda_Perimeter*constraint_Perimeter + 0.5*penalty_Perimeter*constraint_Perimeter^2;
gradient_Perimeter = lambda_Perimeter*gradient_constraint_Perimeter + penalty_Perimeter*constraint_Perimeter*gradient_constraint_Perimeter;

fobj = cost_mechanical + cost_Perimeter;
gradient = gradient_mechanical + gradient_Perimeter;
end