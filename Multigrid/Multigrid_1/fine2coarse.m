function [VCoarseMat] = fine2coarse(fine_nx, fine_ny, VFineMat)
coarse_nx = (fine_nx+1)/2;
coarse_ny = (fine_ny+1)/2;

vCoarseMat = zeros(coarse_nx,coarse_ny);

for i=1:coarse_nx
    for j=1:coarse_ny
        if (i==coarse_nx || i==coarse_nx || j==coarse_ny || j==coarse_ny)
            continue
        else
            VCoarseMat(i,j) = (1.0/16.0)*(VFineMat(2*i-2,2*j) + ...
                VFineMat(2*i,2*j) + VFineMat(2*i-2,2*j-2) + ...
                VFineMat(2*i,2*j-2)) + (1.0/8.0)*(VFineMat(2*i-1,2*j-2) ...
                + VFineMat(2*i,2*j-1) + VFineMat(2*i-1,2*j) + ...
                VFineMat(2*i-2,2*j-1)) + 0.25*VFineMat(2*i-1,2*j-1);
        end
    end
end

end