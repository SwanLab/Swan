function [u_h, gs_iter] = multigrid(nx_coarse, ny_coarse, A_h, b)

initial_guess = zeros(size(b));
conv = 1;
gs_iter = 0;

I_2h_h = construct_I(nx_coarse, ny_coarse);
R_h_2h = construct_R(I_2h_h);
A_2h = R_h_2h * A_h * I_2h_h;

% for eq_num=1:nj_coarse
%     if (id_coarse(eq_num) ~= 0)
%         A_2h(:,eq_num) = 0.0;
%         A_2h(eq_num,:) = 0.0;
%         A_2h(eq_num,eq_num) = 1.0;
%     end
% end

while(conv == 1)
    
    % Apply 3 steps of CG
    %[u_h, conv] = cg_fixed(A_h, b, 3, initial_guess);
    [u_h, conv, local_iter] = gauss_seidel_fixed(A_h, b, initial_guess, 3);
    gs_iter = gs_iter + local_iter;
 
    rh = b - A_h*u_h;
    if (conv == 1)
        disp('Gauss seidel converged')
        return;
    else
        disp('Running multigrid')
        r_2h = R_h_2h * rh;
%         for eq_num=1:nx_coarse*ny_coarse
%             if (id_coarse(eq_num) ~= 0)
%                 r_2h(eq_num) = 0.0;
%             end
%         end

        E_2h = A_2h\r_2h;
        
        E_h = I_2h_h*E_2h;
        
        initial_guess = u_h + E_h;
        conv = 1;
    end
end
end