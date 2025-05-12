classdef VideoGenerator

    properties (Access = private)

        frameRate = 30;
        repeatThreshold = 130;
        moduloSkip = 6;

    end

    methods (Access = public, Static)

        function compute(obj, fileName, totalFrames, dataMatrix, plotFunction)
            v = VideoWriter(fileName);
            v.FrameRate = obj.frameRate;
            open(v);

            figure;
            for i = 1:totalFrames
                s = dataMatrix(i, :);
                plotFunction(s, i);

                frame = getframe(gcf);

                if i <= obj.repeatThreshold
                    repeat = 2;
                elseif mod(i, obj.moduloSkip) ~= 0
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
