classdef VideoGenerator

    methods (Access = public, Static)

        function compute(fileName, totalFrames, repeatThreshold, dataMatrix, plotFunction)

            frameRate = 30;

            v = VideoWriter(fileName);
            v.FrameRate = frameRate;
            open(v);

            figure;
            for i = 1:totalFrames

                if iscell(dataMatrix)
                    s = dataMatrix{i};
                else
                    s = dataMatrix(i, :);
                end

                plotFunction(s, i);

                frame = getframe(gcf);

                if i <= repeatThreshold
                    repeat = 2;
                else
                    repeat = 1;
                end

                for r = 1:repeat
                    writeVideo(v, frame);
                end
                close all
            end

            close(v);
        end
    end
end
