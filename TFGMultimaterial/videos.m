function videos(FRAMES,filename,FrameRate)
%Video Writer
FRAMES(1).cdata = [];
v = VideoWriter(filename);
if isempty(FrameRate)
    v.FrameRate = 30;
else
    v.FrameRate = FrameRate;
end
open(v)
for i=2:length(FRAMES)
    writeVideo(v,FRAMES(i));
end
close(v)