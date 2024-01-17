

function res = residual(u, f, h)
% Calculate the residual of the 2D Poisson equation with Dirichlet BCs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUTS:
% u: 2D array of floats representing the current approximation of u
% f: 2D array of floats representing the right- hand side function f
% h: float representing the grid spacing (dx = dy = h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUTS:
% res: 2D array of floats representing the residual
    
    [m,n ] = size(u);
    res = zeros(m,n);
    
    for i = 2:m-1
        for j = 2:n-1
            res(i,j) = f(i,j) - (1/h^2)*(-4*u(i,j)+ u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1));
        end
    end
    
    %Impose BCs
%     res(1,:) = 0;
%     res(m,:) = 0;
%     res(:,1) = 0;
%     res(:,n) = 0;
%disp(res)
end