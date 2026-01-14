function [on_segment_idx, distances] = points_on_segment(Y, x1, x2)
    % Vector from x1 to x2
    v = x2 - x1;
    v_norm_sq = sum(v.^2);  % Squared length of the segment

    % Vector from x1 to all points in Y
    w = Y - x1;  % Each row is a vector from x1 to yi

    % Project w onto v: compute scalar projection t = dot(w, v) / dot(v, v)
    t = (w * v') / v_norm_sq;  % Column vector of scalars

    % Check if projections fall between 0 and 1 (i.e., within the segment)
    on_segment_mask = (t >= 0) & (t <= 1);

    % Reconstruct projected points
    proj = x1 + t .* v;

    % Check which points are *close enough* to the projected point
    tol = 1e-4;  % Tolerance to account for floating point errors
    distances_to_line = sqrt(sum((Y - proj).^2, 2));
    close_enough_mask = distances_to_line < tol;

    % Final mask: on the segment and close enough to the line
    final_mask = on_segment_mask & close_enough_mask;

    % Return indices and distances
    on_segment_idx = find(final_mask);
    distances = t(final_mask) * sqrt(v_norm_sq);  % distance from x1 along the segment
end
