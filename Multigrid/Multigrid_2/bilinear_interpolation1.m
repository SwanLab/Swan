function prolonged_u = bilinear_interpolation1(u)
% Interpolate a 2D array of grid values using bilinear interpolation

% Get the size of the input array
[m, n] = size(u);

% Create the output array with double the resolution in each dimension
prolonged_u = zeros(2*m, 2*n);

% Copy the original values to the appropriate locations in the output array
prolonged_u(1:2:end, 1:2:end) = u;

% Interpolate in the x-direction
for i = 1:m
    for j = 1:n-1
        prolonged_u(2*i-1, 2*j) = 0.5*(u(i,j) + u(i,j+1));
    end
end

% Interpolate in the y-direction
for i = 1:m-1
    for j = 1:2*n-1
        prolonged_u(2*i, 2*j-1) = 0.5*(u(i,j) + u(i+1,j));
    end
end

% Interpolate in both directions at once
for i = 1:m-1
    for j = 1:n-1
        prolonged_u(2*i, 2*j) = 0.25*(u(i,j) + u(i+1,j) + u(i,j+1) + u(i+1,j+1));
    end
end

end