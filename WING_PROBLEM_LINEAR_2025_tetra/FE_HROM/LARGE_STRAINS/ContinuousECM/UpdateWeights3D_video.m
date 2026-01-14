function UpdateWeights3D_video(DATA,fhandle)

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

%if DATA.ChooseAXES == 1
DATA= DefaultField(DATA,'anglesVIEW',[]) ; 

if isempty(DATA.anglesVIEW)
ok = menu('Set manually the viewpoint for the video; then press ok','ok') ; 
else
    view(DATA.anglesVIEW(1),DATA.anglesVIEW(2)) ; 
end
%if ok == 1
%end

%end

for iter = 1:length(DATA.xALL)
    
    
    UpdateWeights3Dloc(DATA,iter) ;
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if iter == 1
        if  isempty(DATA.anglesVIEW)
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
            'DelayTime',DATA.DelayBetweenFrames);
        else
              imwrite(imind,cm,filename,'gif','WriteMode','append',...
            'DelayTime',6);
        end
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
            'DelayTime',DATA.DelayBetweenFrames);
    end
    
    
end

RotationView1 = linspace(0,360,100) ; 
RotationView2 = linspace(0,0,100) ; 

[ang1,ang2] = view ; 
for irot = 1:length(RotationView1)
    view(ang1+RotationView1(irot),ang2+RotationView2(irot)) ; 
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
     imwrite(imind,cm,filename,'gif','WriteMode','append',...
            'DelayTime',DATA.DelayBetweenFrames);
end



 


% iter_before = max(1,iter-1) ;
% OldWeight = DATA.wALL(:,iter_before) ;
% %
% set(DATA.h_old,'Ydata',OldWeight) ;


% DATA.stepplot = round(Value);
% [DATA]= May29_step(DATA);
% guidata(gcf,DATA) ;
