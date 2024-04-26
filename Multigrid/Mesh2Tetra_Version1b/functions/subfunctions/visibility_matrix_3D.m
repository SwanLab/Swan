function VS=visibility_matrix_3D(V,F)
V=V+rand(size(V))*1e-15;
VS=true([size(V,1) size(V,1)]);
for i=1:size(V,1)
    for j=(i+1):size(V,1)
        Visible=CheckVisiblePoint3D(V,F,i,j);
        VS(i,j)=Visible; VS(j,i)=Visible;
    end
end

