

function u_fine = bilinear_interpolation(u_coarse, u_fine)
% Interpolate solution from coarse to fine grid using bilinear interpolation.

% Extract dimensions of coarse grid
[m_coarse, n_coarse] = size((u_coarse));

% Define dimensions of fine grid
m_fine = 2*m_coarse - 1;
%disp(m_coarse);
%disp(n_coarse);
n_fine = 2*n_coarse - 1;
%disp(n_fine)

% Initialize fine solution
if nargin == 1
    u_fine = zeros(m_fine, n_fine);
end

% Interpolate values
for i = 1:m_coarse-2
    for j = 1:n_coarse-2
        % Coarse indices
        ic = i;
        jc = j;
       
        % Fine indices
        if mod(i, 2) == 0
            i1 = 2*i-1;
            i2 = 2*i+1;
            ic = i/2;
        else
            i1 = 2*i;
            i2 = 2*i;
        end
       
        if mod(j, 2) == 0
            j1 = 2*j-1;
            j2 = 2*j+1;
            jc = j/2;
        else
            j1 = 2*j;
            j2 = 2*j;
            
            
        end
        %disp(j1);
        %disp(j2);
       
        % Interpolate
        if i > 1 && i < m_coarse && j > 1 && j < n_coarse
            u_fine(i1,j1) = u_coarse(i,j);
            u_fine(i2,j1) = 0.5*(u_coarse(i,j) + u_coarse(i+1,j));
            u_fine(i1,j2) = 0.5*(u_coarse(i,j) + u_coarse(i,j+1));
            u_fine(i2,j2) = 0.25*(u_coarse(i,j) + u_coarse(i+1,j) + u_coarse(i,j+1) + u_coarse(i+1,j+1));
            %disp(u_fine);
        else
            u_fine(i1,j1) = u_coarse(i,j);
            u_fine(i2,j1) = u_coarse(i+1,j);
            u_fine(i1,j2) = u_coarse(i,j+1);
            %u_fine(i2,j2) = u_coarse(i+1,j+1);

            if i == 1 
                u_fine(i1,j1) = u_coarse(i,j);
                u_fine(i2,j1) = u_coarse(i+1,j);
            end
             if i == m_coarse 
                u_fine(i1,j1) = u_coarse(i,j);
                u_fine(i2,j1) = u_coarse(i+1,j);
             end
             if j == 1 
                u_fine(i1,j1) = u_coarse(i,j);
                u_fine(i1,j2) = u_coarse(i,j+1);
             end
             if j == n_coarse 
                u_fine(i1,j1) = u_coarse(i,j);
                u_fine(i1,j2) = u_coarse(i,j+1);
            end

        end
        
    end
end
disp(u_fine);

end
