function compact_stress = compact_quadratic_stress(stress)

ngauss = size(stress,1);
nelem = size(stress,3);
compact_stress = zeros(ngauss,6,nelem);
for igauss = 1:ngauss
compact_stress(igauss,1,:) = stress(igauss,1,:).*stress(igauss,1,:);
compact_stress(igauss,2,:) = stress(igauss,2,:).*stress(igauss,2,:);
compact_stress(igauss,3,:) = stress(igauss,3,:).*stress(igauss,3,:);
compact_stress(igauss,4,:) = stress(igauss,2,:).*stress(igauss,3,:);
compact_stress(igauss,5,:) = stress(igauss,1,:).*stress(igauss,3,:);
compact_stress(igauss,6,:) = stress(igauss,1,:).*stress(igauss,2,:);
end



end