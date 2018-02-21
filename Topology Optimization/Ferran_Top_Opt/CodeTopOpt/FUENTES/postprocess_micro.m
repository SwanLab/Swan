function postprocess_micro(Path,newmicros,iter)
    if iter == 1
    %unix(['rm ',Path,'*.res']);
    %unix(['rm ',Path,'*.msh']);
    unix(['rm ',Path,'*.vv']);
    end
    MultipleFiles = '';
    %Path_micros = '/media/aferrer/Alex/Vademecum_20_per1e-3/';
    Path_micros = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/Cases/';
    %newmicros = unique(num2str(floor(newmicros/4.0000001)+1),'rows');
for imicro = 1:length(newmicros)
    
    if newmicros(imicro) <= 2145*4
        micro_ref = floor(newmicros(imicro)/4.0000001)+1; %Micro del primer quadrante
        ncuadrante = newmicros(imicro) - 4*(micro_ref-1);
        %[~,micro_file_msh] = unix(['ls ',Path_micros,num2str(micro_ref),'/*.msh']);
        %[~,micro_file_res] = unix(['ls ',Path_micros,num2str(micro_ref),'/*.res']);
        [~,micro_file_msh] = unix(['find ',Path_micros,num2str(micro_ref),'/*.msh -not -iname "ElasticMicro_0.flavia.msh"']);
        [~,micro_file_res] = unix(['find ',Path_micros,num2str(micro_ref),'/*.res -not -iname "ElasticMicro_0.flavia.res"']);
        unix(['cp ',micro_file_msh(1:end-1),[' ',Path,num2str((imicro)),'_',num2str(iter),'.msh']]);
        unix(['cp ',micro_file_res(1:end-1),[' ',Path,num2str((imicro)),'_',num2str(iter),'.res']]);
        changeTimeStepNum(str2double([num2str(imicro),num2str(iter)]),[[Path,num2str((imicro))],'_',num2str(iter),'.res'])
        MultipleFiles = [MultipleFiles,[Path,num2str((imicro))],'_',num2str(iter),'.res '];
        
        if ncuadrante ~= 1 % Rotation o mirroring of micro del mfl
            change_msh_file([Path,num2str((imicro)),'_',num2str(iter),'.msh'],ncuadrante)
        end
        
    else %Tomar circulo como micro
        unix(['cp /home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/circuloV0_6.msh',[' ',Path,num2str((imicro)),'_',num2str(iter),'.msh']]);
        unix(['cp /home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/circuloV0_6.res',[' ',Path,num2str((imicro)),'_',num2str(iter),'.res']]);
        changeTimeStepNum(str2double([num2str(imicro),num2str(iter)]),[[Path,num2str((imicro))],'_',num2str(iter),'.res'])
        MultipleFiles = [MultipleFiles,[Path,num2str((imicro))],'_',num2str(iter),'.res '];
    end
    
    
    
    
end
writebash(Path,MultipleFiles)
%opengid
end

function changeTimeStepNum(imicro,Path)
fid  = fopen(Path);
f=fread(fid,'*char')';
fclose(fid);
timestep = f(282:286);
f = strrep(f,timestep,['     ',num2str(imicro)]);

fid  = fopen(Path,'w');
fprintf(fid,'%s',f);
fclose(fid);
end

function change_msh_file(Path,ncuadrante)
fid = fopen(Path,'r');
C = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
C = C{1};
%C = textread([Path,'RVE04N3_11.flavia.msh'], '%s','delimiter', '\n');

Cnew = C;

lines = [7:3287];
coor = zeros(length(lines),3);
newcoor = zeros(size(coor));
for iline = 1:length(lines)
    
    coor(iline,:) = str2num(C{lines(iline)});
    newcoor(iline,1) = coor(iline,1);
    switch ncuadrante
        case 2  %% pi/4 < phi < 3pi/4
            P = [0 -1;-1 0]; t = [1 1]';
        case 3  %% 3pi/4 < phi < 5pi/4
            P = [-1 0; 0 1]; t = [1 0]';
        case 4  %% 5pi/4 < phi < 7pi/4
            P = [0 1; -1 0]; t = [0 1]';
    end
    newcoor(iline,2:3) = P*coor(iline,2:3)' + t;
    Cnew{lines(iline)} = num2str(newcoor(iline,:));

end

fid =  fopen(Path,'w+');
fprintf(fid, '%s\n', Cnew{:});
fclose(fid);

end

function writebash(Path,MultipleFiles)
fID = fopen([Path,'/my_batch.bch'],'w+');
fprintf(fID,['Postprocess \n']);
fprintf(fID,['Files ReadMultiple ','{',MultipleFiles,'}',' \n']);
fprintf(fID,[char(39),'Zoom Frame']);
fclose(fID);
end