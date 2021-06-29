function squeezed = squeezeParticular(t,dim)
z = size(t);
index = setdiff(1:length(z),dim);
squeezed = reshape(t,[z(index) z(dim)]);
end