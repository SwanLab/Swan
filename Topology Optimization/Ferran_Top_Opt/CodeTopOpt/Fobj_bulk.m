function F_obj_and_theta_vademecum
Path = '/media/Alex/Vademecum';
fid = fopen([Path,'/',num2str(isnap),'/executed_cases.txt'],'r');
C = textscan(lastline,repmat(['%s%f'],1,7));
fclose(fid);  
end

%
%C = textscan(lastline, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
%    k = 1;
%    for i = 2:2:28
%        Datos{k} = [C{i}];
%        k = k + 1;
%    end
%
%
%figure(1)
%  plot(A);
%figure(2)
%plot(theta);
%
%
%end
%
