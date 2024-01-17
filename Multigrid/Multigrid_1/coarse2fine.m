function [VFineMat] = coarse2fine(coarse_nx, coarse_ny, VCoarseMat)

fine_nx = (coarse_nx+1)/2;
fine_ny = (coarse_ny+1)/2;
VFineMat = reshape(fine_nx, fine_ny);

for i=0:coarse_nx
    for j=0:coarse_ny
         VFineMat(2*i+1,2*j+1) = VCoarseMat(i+1,j+1);
    end
end
for i=0:coarse_nx
    for j=0:coarse_ny-1
         VFineMat(2*i+1, 2*j+2) = 0.5*(VCoarseMat(i+1,j+1) + ...
             VCoarseMat(i+1,j+2));
    end
end
for i=0:coarse_nx-1
    for j=0:coarse_ny
         VFineMat(2*i+2,2*j+1) = 0.5*(VCoarseMat(i+1,j+1) + ...
             VCoarseMat(i+2,j+1));
    end
end
for i=0:coarse_nx-1
    for j=0:coarse_ny-1
         VFineMat(2*i+2,2*j+2) = 0.25*(VCoarseMat(i+1,j+1) ...
             + VCoarseMat(i+2,j+1) + VCoarseMat(i+1,j+2) + ...
             VCoarseMat(i+2,j+2));
    end
end

end

