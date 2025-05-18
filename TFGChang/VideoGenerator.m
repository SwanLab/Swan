classdef VideoGenerator

    methods (Access = public, Static)

        function compute(fileName, totalFrames, dataMatrix, plotFunction)
            
            frameRate = 30;
            repeatThreshold = 30;
            moduloSkip = 6;
            
            v = VideoWriter(fileName);
            v.FrameRate = frameRate;
            open(v);

            figure;
            for i = 1:totalFrames
                s = dataMatrix(i, :);
                plotFunction(s, i);

                frame = getframe(gcf);

                if i <= repeatThreshold
                    repeat = 6;
                elseif mod(i, moduloSkip) ~= 0
                    continue;
                else
                    repeat = 1;
                end

                for r = 1:repeat
                    writeVideo(v, frame);
                end
            end

            close(v);
        end
    end
end
