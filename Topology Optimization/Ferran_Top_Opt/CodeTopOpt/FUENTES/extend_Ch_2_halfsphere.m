function [C_total,def_total,phi_total,theta_total,phi_matrix_total_unique,theta_matrix_total,number_matrix_total_unique] = extend_Ch_2_halfsphere(Ch_all,phi_all,theta_all,phi_matrix,theta_matrix,number_matrix)
n_cuadrantes = 4;
nsnap = size(Ch_all,1);
C_total = zeros(n_cuadrantes*nsnap,size(Ch_all,2));
def_total = zeros(n_cuadrantes*nsnap,3);
phi_total = zeros(n_cuadrantes*nsnap,1);
theta_total = zeros(n_cuadrantes*nsnap,1);

phi_matrix_total = zeros(size(phi_matrix,1)*n_cuadrantes,size(phi_matrix,1));
number_matrix_total = zeros(size(number_matrix,1)*n_cuadrantes,size(number_matrix,1));
theta_matrix_total = theta_matrix;

for isnap = 1:nsnap 
    phi = phi_all(isnap);
    theta = theta_all(isnap);
    Ch = Ch_all(isnap,:);
    Ch_mat = [Ch(1) Ch(2) Ch(3); Ch(2) Ch(4) Ch(5); Ch(3) Ch(5) Ch(6)];
    
   % if phi >= 0
        
%         phi_all4 = [phi, pi/2 - phi, phi, pi/2 - phi]';
%         theta_all4 = [theta, theta, pi-theta, pi-theta]';
%         def_all4 = phitheta2strain(phi_all4,theta_all4);
        

         phi_all4 = [phi,pi/2 - phi,pi + phi, -pi/2 - phi]';
         theta_all4 = [theta, theta, theta, theta]';
         def_all4 = phitheta2strain(phi_all4,theta_all4);

         
         [row,column] = find(number_matrix == isnap);
         rows_total = column*(0:n_cuadrantes-1)+ row;
         rows_total = [row,2*column - row+1, 2*column + row, 4*column - row+1];
         phi_matrix_total(rows_total,column) = phi_all4;
         number_matrix_total(rows_total,column) = n_cuadrantes*(isnap-1)+1:n_cuadrantes*isnap; 
         %          figure(1)
%          pause(1)
%          spy(phi_matrix_total)
%         T(:,:,1) = eye(3);
%         T(:,:,2) = [0 sign(phi) 0; sign(phi) 0 0; 0 0 1]; 
%         T(:,:,3) = [1 0 0; 0 1 0; 0 0 -1];
%         T(:,:,4) = [0 sign(phi) 0; sign(phi) 0 0; 0 0 -1]; 

        T(:,:,1) = eye(3);
        T(:,:,2) = [0 1 0; 1 0 0; 0 0 1]; 
        T(:,:,3) = [1 0 0; 0 1 0; 0 0 -1];
        T(:,:,4) = [0 -1 0; -1 0 0; 0 0 1]; 



        
        C_all4 = zeros(n_cuadrantes,size(Ch_all,2));
        for i =1:n_cuadrantes
            Cquart = T(:,:,i)'*Ch_mat*T(:,:,i);
            C_all4(i,:) = [Cquart(1,1) Cquart(1,2) Cquart(1,3) Cquart(2,2) Cquart(2,3) Cquart(3,3)]; 
        end
        C_total(n_cuadrantes*(isnap-1)+(1:n_cuadrantes),:) = C_all4;
        def_total(n_cuadrantes*(isnap-1)+(1:n_cuadrantes),:) = def_all4;
        phi_total(n_cuadrantes*(isnap-1)+(1:n_cuadrantes),1) = phi_all4;
        theta_total(n_cuadrantes*(isnap-1)+(1:n_cuadrantes),1) = theta_all4;
  %  end

    
end
%C_total = unique(C_total,'rows');
phi_matrix_total_unique = zeros(size(phi_matrix_total));
number_matrix_total_unique = zeros(size(number_matrix_total));
ncolumn = length(theta_matrix_total); 
for icolumn = 1:ncolumn; 
    c = mod(phi_matrix_total(1:icolumn*n_cuadrantes,icolumn),2*pi);
    [a,b] = sort(c); 
    phi_matrix_total(1:icolumn*n_cuadrantes,icolumn) = a;
    number_matrix_total(1:icolumn*n_cuadrantes,icolumn) = number_matrix_total(b,icolumn);
    c = phi_matrix_total(1:icolumn*n_cuadrantes,icolumn);
    num_col_unique = number_matrix_total(1:icolumn*n_cuadrantes,icolumn);
    

    
    [a,b] = uniquetol(c,1e-4);
    
    %[a,b] = unique(c);
    %b = ~any(triu(abs(bsxfun(@minus,c,c.'))<1e-3,1));
    %a = c(b);
    phi_matrix_total_unique(1:length(a),icolumn) = a;
    number_matrix_total_unique(1:length(a),icolumn) = num_col_unique(b);
    
    
    
end





end