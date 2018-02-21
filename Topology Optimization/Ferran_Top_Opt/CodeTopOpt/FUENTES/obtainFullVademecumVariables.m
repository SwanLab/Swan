function [cuadrante,coord_full,Ch_full,ener_full,phi_full,theta_full,strain_full] = obtainFullVademecumVariables(coordgid,Ch,ener,phi,theta,strain)
% Compute information of the whole sphere
nSnap = size(coordgid,1);
color = ['bkrgycbkr'];
rounding = 5;

%rotation matrix on the sphere
cuadrante(1).R = [1 0 0; 0 1 0;0 0 1];
cuadrante(2).R = [0 1 0; 1 0 0; 0 0 1];
cuadrante(3).R = [-1 0 0; 0 -1 0; 0 0 1];
cuadrante(4).R = [0 -1 0; -1 0 0; 0 0 1];
cuadrante(5).R = [1 0 0; 0 1 0;0 0 -1];
cuadrante(6).R = [0 1 0; 1 0 0; 0 0 -1];
cuadrante(7).R = [-1 0 0; 0 -1 0; 0 0 -1];
cuadrante(8).R = [0 -1 0; -1 0 0; 0 0 -1];

%tranformation of the topology
cuadrante(1).P = [1 0 ; 0 1];
cuadrante(2).P = [0 -1;-1 0];
cuadrante(3).P = [-1 0; 0 1];
cuadrante(4).P = [0 1; -1 0];
cuadrante(5).P = cuadrante(3).P;
cuadrante(6).P = cuadrante(4).P;
cuadrante(7).P = cuadrante(1).P;
cuadrante(8).P = cuadrante(2).P;

cuadrante(1).t = [0 0]';
cuadrante(2).t = [1 1]';
cuadrante(3).t = [1 0]';
cuadrante(4).t = [0 1]';
cuadrante(5).t = cuadrante(3).t;
cuadrante(6).t = cuadrante(4).t;
cuadrante(7).t = cuadrante(1).t;
cuadrante(8).t = cuadrante(2).t;

cuadrante(1).color = color(1);
cuadrante(2).color = color(2);
cuadrante(3).color = color(3);
cuadrante(4).color = color(4);
cuadrante(5).color = color(5);
cuadrante(6).color = color(6);
cuadrante(7).color = color(7);
cuadrante(8).color = color(8);

cuadrante(1).a_phi = 1;
cuadrante(2).a_phi = -1;
cuadrante(3).a_phi = 1;
cuadrante(4).a_phi = -1;
cuadrante(5).a_phi = cuadrante(1).a_phi;
cuadrante(6).a_phi = cuadrante(2).a_phi;
cuadrante(7).a_phi = cuadrante(3).a_phi;
cuadrante(8).a_phi = cuadrante(4).a_phi;

cuadrante(1).b_phi = 0;
cuadrante(2).b_phi = pi/2;
cuadrante(3).b_phi = pi;
cuadrante(4).b_phi = 3*pi/2;
cuadrante(5).b_phi = cuadrante(1).b_phi;
cuadrante(6).b_phi = cuadrante(2).b_phi;
cuadrante(7).b_phi = cuadrante(3).b_phi;
cuadrante(8).b_phi = cuadrante(4).b_phi;

cuadrante(1).a_theta = 1;
cuadrante(2).a_theta = 1;
cuadrante(3).a_theta = 1;
cuadrante(4).a_theta = 1;
cuadrante(5).a_theta = -1;
cuadrante(6).a_theta = -1;
cuadrante(7).a_theta = -1;
cuadrante(8).a_theta = -1;

cuadrante(1).b_theta = 0;
cuadrante(2).b_theta = 0;
cuadrante(3).b_theta = 0;
cuadrante(4).b_theta = 0;
cuadrante(5).b_theta = cuadrante(1).b_theta;
cuadrante(6).b_theta = cuadrante(2).b_theta;
cuadrante(7).b_theta = cuadrante(3).b_theta;
cuadrante(8).b_theta = cuadrante(4).b_theta;



cuadrante(1).coord = round(coordgid(:,2:end),rounding);
cuadrante(1).strain = strain;
cuadrante(1).Ch = Ch;
cuadrante(1).ener = ener;
cuadrante(1).phi = phi;
cuadrante(1).theta = theta;
cuadrante(1).pointer_2_first_cuadrante = 1:size(coordgid,1); %pointer al primer cuadrante
cuadrante(1).num_global = 1:size(coordgid,1); %numeracion global
cuadrante(1).npoint = size(coordgid,1);


coord_full = cuadrante(1).coord;
strain_full = cuadrante(1).strain;
Ch_full = cuadrante(1).Ch;
ener_full = cuadrante(1).ener;
phi_full = cuadrante(1).phi;
theta_full = cuadrante(1).theta;
pointer_2_first_cuadrante_full = cuadrante(1).pointer_2_first_cuadrante;

%plot3(coord_full(:,1),coord_full(:,2),coord_full(:,3),'b+');
%hold on


for icuadrante = 2:numel(cuadrante)
    coord_cuadrante = [];
    for ipoint = 1:nSnap
            coord_cuadrante(ipoint,:) = (cuadrante(icuadrante).R*cuadrante(1).coord(ipoint,:)')';
            strain_cuadrante(ipoint,:) = (cuadrante(icuadrante).R*cuadrante(1).strain(ipoint,:)')';
            Ch_cuadrante(:,:,ipoint) = cuadrante(icuadrante).R'*cuadrante(1).Ch(:,:,ipoint)*cuadrante(icuadrante).R;
            ener_cuadrante(ipoint) = cuadrante(1).ener(1,ipoint);
            phi_cuadrante(ipoint) = cuadrante(icuadrante).a_phi*cuadrante(1).phi(1,ipoint) + cuadrante(icuadrante).b_phi;
            theta_cuadrante(ipoint) = cuadrante(icuadrante).a_theta*cuadrante(1).theta(1,ipoint) + cuadrante(icuadrante).b_theta;
            pointer_2_first_cuadrante_cuadrante(ipoint) = ipoint;
           % ener_cuadrante(ipoint) = 
    end
    coord_full_repeted_augmented = cat(1,coord_full,coord_cuadrante);
    strain_full_repeted_augmented = cat(1,strain_full,strain_cuadrante);
    Ch_full_repeted_augmented = cat(3,Ch_full, Ch_cuadrante);
    ener_full_repeted_augmented = cat(2,ener_full,ener_cuadrante);
    phi_full_repeted_augmented = cat(2,phi_full,phi_cuadrante);
    theta_full_repeted_augmented = cat(2,theta_full,theta_cuadrante);
    pointer_2_first_cuadrante_full_repeted_augmented = cat(2,pointer_2_first_cuadrante_full,pointer_2_first_cuadrante_cuadrante);
    [coord_cuadrante_unique,index_coord] = setdiff(coord_full_repeted_augmented,coord_full,'rows');
    
    cuadrante(icuadrante).coord = coord_cuadrante_unique;
    cuadrante(icuadrante).strain = strain_full_repeted_augmented(index_coord,:);
    cuadrante(icuadrante).Ch = Ch_full_repeted_augmented(:,:,index_coord);
    cuadrante(icuadrante).ener = ener_full_repeted_augmented(1,index_coord);
    cuadrante(icuadrante).phi = phi_full_repeted_augmented(1,index_coord);
    cuadrante(icuadrante).theta = theta_full_repeted_augmented(1,index_coord);
    cuadrante(icuadrante).pointer_2_first_cuadrante = pointer_2_first_cuadrante_full_repeted_augmented(1,index_coord);
    cuadrante(icuadrante).npoint = size(coord_cuadrante_unique,1);
    cuadrante(icuadrante).num_global = [1:cuadrante(icuadrante).npoint] + cuadrante(icuadrante-1).num_global(end);
    
   % plot3(coord_cuadrante_unique(:,1),coord_cuadrante_unique(:,2),coord_cuadrante_unique(:,3),[color(icuadrante),'+'])
    coord_full = cat(1,coord_full, coord_cuadrante_unique);
    strain_full = cat(1,strain_full, strain_full_repeted_augmented(index_coord,:));
    Ch_full = cat(3,Ch_full, Ch_full_repeted_augmented(:,:,index_coord));
    ener_full = cat(2,ener_full,ener_full_repeted_augmented(1,index_coord));
    phi_full = cat(2,phi_full,phi_full_repeted_augmented(1,index_coord));
    theta_full = cat(2,theta_full,theta_full_repeted_augmented(1,index_coord));
    pointer_2_first_cuadrante_full = cat(2,pointer_2_first_cuadrante_full,pointer_2_first_cuadrante_full_repeted_augmented(1,index_coord));
end

end