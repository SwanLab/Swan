function [fobj,gradient] = assamble_cost_gradient_mechanical(compliance,gradient_compliance,constraint_volume,gradient_constraint_volume,...
                                                            lambda_volume,penalty_volume)

cost_volume = lambda_volume*constraint_volume + 0.5*penalty_volume*constraint_volume^2;

gradient_volume = lambda_volume*gradient_constraint_volume + penalty_volume*constraint_volume*gradient_constraint_volume;

fobj = compliance + cost_volume;
gradient = gradient_compliance + gradient_volume;
end