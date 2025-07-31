function makePolicyGif(weight_history, params, getActiveTiles, filename)
    if isfile(filename)
        delete(filename);
    end

    for i = 1:size(weight_history, 2)
        ep = i * params.gifFreq;
        weights = weight_history(:, i);

        img = gif_plotPolicyColor2D(weights, params, getActiveTiles, ep);
        [A, map] = rgb2ind(img, 256);

        if i == 1
            imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
        else
            imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
        end
    end
end