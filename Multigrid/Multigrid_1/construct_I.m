function [global_I] = construct_I(nx_coarse, ny_coarse)

mg_weight_mat = ...
    [1 0 0 0;
    0.5 0.5 0 0;
    0 1 0 0;
    0.5 0 0.5 0;
    0.25 0.25 0.25 0.25;
    0 0.5 0 0.5;
    0 0 1 0;
    0 0 0.5 0.5;
    0 0 0 1];

nx_fine = 2*nx_coarse -1;
ny_fine = 2*ny_coarse -1;

global_I = zeros(nx_fine*ny_fine, nx_coarse*ny_coarse);
for j = 1:ny_coarse-1
    for i = 1:nx_coarse-1 
        
        % Obtain global ids of coarse mesh of current cell
        cbl = (j-1)*nx_coarse + i;
        coarse_indices = [cbl, cbl+1, cbl+nx_coarse, cbl+nx_coarse+1];
        
        fbl = (2*j-2)*nx_fine + 2*i - 2;
        for k=1:3
            for l=1:3
                % Get global node number in finer mesh
                % This will be the row number of global_I
                global_I(fbl + (k-1)*nx_fine + l, coarse_indices) = ...
                    global_I(fbl + (k-1)*nx_fine + l, coarse_indices) + ...
                    mg_weight_mat((k-1)*3 + l, :);        
            end
        end
    end
end
global_I = global_I./sum(global_I,2);

% R = zeros(nx_coarse*ny_coarse,nx_fine*ny_fine);
% for j = 1:2:ny_fine-1
%     for i = 1:2:nx_fine-1 
%         
%         % Obtain global ids of coarse mesh of current cell
%         fbl = (j-1)*nx_fine + i;
%         fine_indices = [fbl, fbl+1, fbl+2, ...
%                         fbl+nx_fine, fbl+nx_fine+1, fbl+nx_fine+2,...
%                         fbl+2*nx_fine, fbl+2*nx_fine+1, fbl+2*nx_fine+2];
%         
%         cbl = (j-1)*nx_coarse/2 + (i-1)/2;
%         for k=1:2
%             for l=1:2
%                 % Get global node number in finer mesh
%                 % This will be the row number of global_I
%                 R(cbl + (k-1)*nx_coarse + l, fine_indices) = ...
%                     R(cbl + (k-1)*nx_coarse + l, fine_indices) + ...
%                     mg_weight_mat2((k-1)*2 + l, :);        
%             end
%         end
%     end
% end
% % R = R./sum(R,2);

% for i = 1:nx_coarse
%     global_I(i,:)=0;
%     global_I(i,i)=1;
%     global_I(nx_coarse*nx_coarse-i-1,:)=0;
%     global_I(nx_coarse*nx_coarse-i-1,nx_coarse*nx_coarse-i-1)=1;
% end
% 
% for i = 1:nx_coarse:nx_coarse*nx_coarse
%     global_I(i,:)=0;
%     global_I(i,i)=1;
%     global_I(i+nx_coarse-1,:)=0;
%     global_I(i + nx_coarse -1,i+nx_coarse-1)=1;
% 
% end
    
end