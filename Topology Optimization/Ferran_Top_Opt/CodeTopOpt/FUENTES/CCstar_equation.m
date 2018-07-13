function [C_Cstar,gradient_constraint] = CCstar_equation (design_variable,compliancefun)

global post_info

% Re-compute equilibrium update and inform user
if ~isequal(post_info.xold,design_variable)
    compliancefun(design_variable);
    message = 'Re-computed equilibrium update due to C-C* constraint';
    try
        cprintf('err','%s\n',message);
    catch ME
        warning('cprintf failed!');
        fprintf('%s\n',message);
    end
end

Ch_star_div = post_info.Ch_star;
Ch_star_div (abs(Ch_star_div) < 1e-5) = 1;
C_C = (post_info.structural_values.matCh - post_info.Ch_star)./Ch_star_div;

% C-C*
sq2 = sqrt(2);
weights = [1,1,1,sq2,sq2,sq2]';
C_Cstar = weights.*[C_C(1,1); ...
                   C_C(2,2); ...
                   C_C(3,3); ...
                   C_C(2,3); ...
                   C_C(1,3); ...
                   C_C(1,2)];
       
problembsc.costfunct = 'linear_components';
neq = 6;
selectiveC_Cstar = zeros(3,3,neq);

% Eqn 1
selectiveC_Cstar(1,1,1) = 1;
selective_Ch_star_div(1) = Ch_star_div(1,1);

% Eqn 2
selectiveC_Cstar(2,2,2) = 1;
selective_Ch_star_div(2) = Ch_star_div(2,2);

% Eqn 3
selectiveC_Cstar(3,3,3) = 1;
selective_Ch_star_div(3) = Ch_star_div(3,3);

% Eqn 4
selectiveC_Cstar(2,3,4) = 1;
selective_Ch_star_div(4) = Ch_star_div(2,3);

% Eqn 5
selectiveC_Cstar(1,3,5) = 1;
selective_Ch_star_div(5) = Ch_star_div(1,3);

% Eqn 6
selectiveC_Cstar(1,2,6) = 1;
selective_Ch_star_div(6) = Ch_star_div(1,2);

gradient_constraint = zeros(post_info.dim.nelem,neq);
for i = 1:neq
    problembsc.selectiveC_Cstar = selectiveC_Cstar(:,:,i)./selective_Ch_star_div(i);
    gradient_constraint(:,i) = weights(i)*compliance_derivative_micro([],post_info.DtC,problembsc,post_info.dim,post_info.ngaus);
end



end