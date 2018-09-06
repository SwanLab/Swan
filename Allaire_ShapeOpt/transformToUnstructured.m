function [vect] = transformToUnstructured(mat)
load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));

if length(mat) == length(A1)
    A = A1;
else
    A = A0;
end

if length(dim) == 2
    vect(A(:,:)) = mat(:,:);
else
    vect(A(:,:,:)) = mat(:,:,:);
end
vect = vect';
end

