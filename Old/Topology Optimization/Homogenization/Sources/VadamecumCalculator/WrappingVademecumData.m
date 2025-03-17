function WrappingVademecumData
% 
% for i = 2:20
%    RunningVademecumInParalel(i); 
% end
path = '/media/alex/My Passport/Vademecum/';
cellVar = cell(20,20,10);
for iMx = 1:20
    fileName = 'OptimalSuperEllipseDifferentStressSign';
    d = load([path,fileName,num2str(iMx),'.mat']);
    for iMy = 1:20
        for iPhi = 1:10
            cellVar{iMx,iMy,iPhi} = d.c{iMx,iMy,iPhi};
        end
    end    
end

end
