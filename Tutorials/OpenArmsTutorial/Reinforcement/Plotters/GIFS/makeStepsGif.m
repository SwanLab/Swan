function makeStepsGif(weight_history, env, params, filename)
    
    for i = 1:size(weight_history, 2)
        episode = i * params.gifFreq;
        w = weight_history(:, i);
        img = gif_plotStepsColor2D(w, env, params, episode);
        [A, map] = rgb2ind(img, 256);

        if i == 1
            imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
        else
            imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
        end
    end
end