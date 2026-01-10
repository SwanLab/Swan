function [K_axial,K_bending_z,MSG] = PropStiffBeam2D(Kskel_11,Kskel_12,Kskel,L,MSG)


   MSG{end+1}='----------------------------';
        %MSG{end+1}='Axial stiffness (E*A)';
        K_axial = Kskel_11(1,1)*L ;
        MSG{end+1}=['K_axial = EA =',num2str(K_axial,5)];
        MSG{end+1}='----------------------------';
        %MSG{end+1}='Bending stiffness -z (E*I_z)';
        nrows = size(Kskel_11,1) +size(Kskel_12,1) ;
        DISP = zeros(nrows,1);
        DISP(3) = 1;
        DISP(6) = -1;
        FORCES = Kskel*DISP ;
        %FORCES(6) Should be equal to 4EI/L-2EI/L  = 2EI/L
        EI = abs(FORCES(3))*L/2 ;  % Notice that this coincide with Kskel(6,6)
        K_bending_z = EI ;
        MSG{end+1}=['K_bending-z = EI_z =',num2str(K_bending_z,5)];
        
%         
%         if  CALCULATE_SHEAR_STIFFNESS == 0
%             MSG{end+1}='----------------------------';
%             MSG{end+1}='Shear stiffness depends on the length of the slice')
%             MSG{end+1}='Here we simply display component K(2,2) and K(3,3) multiplied by the length')
%             MSG{end+1}=['K_SHEAR_y  = ',num2str(Kskel(2,2)*L)])
%             MSG{end+1}=['K_SHEAR_z  = ',num2str(Kskel(3,3)*L)])
%             
%         else
%             K_11_inv = inv(Kskel_11) ;
%             
%             FORCES = zeros(6,1);
%             FORCES(2) = 1;
%             %FORCES(6) = -FORCES(2)*L ;
%             a = K_11_inv*FORCES ;
%             MSG{end+1}='----------------------------')
%             MSG{end+1}='Shear stiffness: only valid for symmetric cross-sections')
%             MSG{end+1}='Besides, it depends on the actual length of the slice')
%             MSG{end+1}='----------------------------')
%             % No we compute the displacement induced considering only the effect of
%             % bending
%             K_11_bending = zeros(6,6) ;
%             K_11_bending(1,1) = K_axial/L ;
%             K_11_bending(4,4) = K_torsion/L ;
%             
%             K_11_bending(2,2) = 12*K_bending_z/L^3 ;
%             K_11_bending(2,6) = 6*K_bending_z/L^2 ;
%             K_11_bending(3,3) = 12*K_bending_y/L^3 ;
%             K_11_bending(3,5) = 6*K_bending_y/L^2 ;
%             
%             K_11_bending(6,6) = 4*K_bending_z/L ;
%             K_11_bending(5,5) = 4*K_bending_y/L ;
%             
%             a_bending = K_11_bending\FORCES ;
%             % Therefore
%             a_shear = a-a_bending ;
%             K_shear = L*FORCES(2)/a_shear(2) ;
%             MSG{end+1}=['K_shear_y = G A_ry =',num2str(K_shear,5)])
%             
%             
%             FORCES = zeros(6,1);
%             FORCES(3) = 1;
%             a = K_11_inv*FORCES ;
%             MSG{end+1}='----------------------------')
%             a_bending = K_11_bending\FORCES ;
%             a_shear = a-a_bending ;
%             
%             K_shear = L*FORCES(3)/a_shear(3) ;
%             
%             
%             MSG{end+1}=['K_shear_z = G A_ry =',num2str(K_shear,5)])
%             
%         end