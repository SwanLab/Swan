function UpdateWeights1Ddecm_video(DATA,fhandle)

if nargin == 0
    load('tmp.mat')
end


% fhandle = guihandles(gcf) ; % --> Tag
% hread = getfield(fhandle,'slide_tag');
% VAL = get(hread,'Value') ;
figure(fhandle)
filename = DATA.NameVideo;
% v = [-5 -2 5];
% [caz,cel] = view(v)

%ok = menu('Set manually the viewpoint for the video; then press ok','ok') ; 

%if ok == 1
%end

for iter = 1:size(DATA.wALL,2)
    
    UpdateWeights1DlocDECM(DATA,iter) ; 

  %  UpdateWeights1Dloc_DECM(DATA,iter) ;
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if iter == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
            'DelayTime',DATA.DelayBetweenFrames);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
            'DelayTime',DATA.DelayBetweenFrames);
    end
    
    
end
