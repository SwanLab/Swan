function value_sequence = sequence_of_stopping_Criteria_parameter(value_first,value_last,epsilon_iter)
ne = length(epsilon_iter);
theta_frac = (value_first/value_last)^(1/(ne-1));
value_sequence = value_first./(theta_frac.^(0:1:(ne-1)));
end