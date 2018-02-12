function [C_Ciso,gradient_constraint] = isotropy_equation2D_scaled (design_variable,compliancefun)

global post_info

% Re-compute equilibrium update and inform user
if ~isequal(post_info.xold,design_variable)
    compliancefun(design_variable);
    message = 'Re-computed equilibrium update due to Isotropy constraint';
    try
        cprintf('err','%s\n',message);
    catch ME
        warning('cprintf failed!');
        fprintf('%s\n',message);
    end
end

matCh = post_info.structural_values.matCh;

% C-Ciso    
C_Ciso = [1 - matCh(2,2)/matCh(1,1); ...
          -1*(1 - matCh(2,2)/matCh(1,1)); ...
          matCh(3,3)/matCh(1,1) - 0.5*(1 + matCh(2,2)/matCh(1,1)) - matCh(2,1)/matCh(1,1); ...
          matCh(2,3); ...
          matCh(1,3); ...
          0];

problembsc.costfunct = 'linear_components';
neq = 6;
selectiveC_Cstar = zeros(3,3,neq);

% Eqn 1
selectiveC_Cstar(1,1,1) = matCh(2,2)/(matCh(1,1)^2);
selectiveC_Cstar(2,2,1) = -1/matCh(1,1);

% Eqn 2
selectiveC_Cstar(1,1,2) = -matCh(2,2)/(matCh(1,1)^2);
selectiveC_Cstar(2,2,2) = 1/matCh(1,1);

% Eqn 3
selectiveC_Cstar(1,1,3) = -matCh(3,3)/(matCh(1,1)^2) + 0.5*matCh(2,2)/(matCh(1,1)^2) + matCh(2,1)/(matCh(1,1)^2);
selectiveC_Cstar(2,2,3) = -0.5/matCh(1,1);
selectiveC_Cstar(3,3,3) = 1/matCh(1,1);
selectiveC_Cstar(2,1,3) = -1/matCh(1,1);

% Eqn 4
selectiveC_Cstar(2,3,4) = 1;

% Eqn 5
selectiveC_Cstar(1,3,5) = 1;

% Eqn 6
% all zeros

gradient_constraint = zeros(post_info.dim.nelem,neq);
for i = 1:neq
    problembsc.selectiveC_Cstar = selectiveC_Cstar(:,:,i);
    gradient_constraint(:,i) = compliance_derivative_micro([],post_info.DtC,problembsc,post_info.dim,post_info.ngaus);
end



end