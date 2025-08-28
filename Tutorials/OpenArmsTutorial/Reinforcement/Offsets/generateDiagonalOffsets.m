function offsets = generateDiagonalOffsets(n_dims, n_tilings)
    base = linspace(0, 1, n_tilings + 1);
    base(end) = [];  % remove the last to avoid overlap
    offsets = repmat(base, n_dims, 1);
end