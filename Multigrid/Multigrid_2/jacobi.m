

function u_new = jacobi(u, f, h, omega, Npre, Npost)
% Performs Npre pre-smoothing and Npost post-smoothing Jacobi iterations on u

% Initialize solution guess
    u_new = u;
    
    % Pre-smoothing
    for i = 1:Npre
        u_old = u_new;
        for j = 2:size(u,1)-1
            for k = 2:size(u,2)-1
                u_new(j,k) = (1-omega)*u_old(j,k) + omega*(0.25*(u_old(j-1,k) + u_old(j+1,k) + u_old(j,k-1) + u_old(j,k+1)) - h^2*f(j,k));
            end
        end
    end
    
    % Post-smoothing
    for i = 1:Npost
        u_old = u_new;
        for j = 2:size(u,1)-1
            for k = 2:size(u,2)-1
                u_new(j,k) = (1-omega)*u_old(j,k) + omega*(0.25*(u_old(j-1,k) + u_old(j+1,k) + u_old(j,k-1) + u_old(j,k+1)) - h^2*f(j,k));
            end
        end
    end

end