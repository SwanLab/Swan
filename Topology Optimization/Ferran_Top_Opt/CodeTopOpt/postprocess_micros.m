function postprocess_micros(Path_micros,micros_evol,case_post,macro_changing)
    
    if macro_changing
        macro_changing_value = 'macro_chang';
    else
        macro_changing_value = 'no_macro_chang';
    end



    for imicro = 1:size(micros_evol,1)
        file_tcl_name = [Path_micros,'tcl_Video_file_micro_',num2str(imicro),'.tcl'];
        write_tcl_file(file_tcl_name,Path_micros,[Path_micros,num2str(imicro)],size(micros_evol,2),case_post,macro_changing_value)
%       unix(['/opt/GiDx64/12.1.2d/gid_offscreen -t "source ',file_tcl_name,'"'])
       % unix(['cp /home/aferrer/Documents/Doctorat/Tesi/FEM_DT/createGidFigures.tcl ',Path_micros])
        if macro_changing
        macro_changing_value = 'macro_chang';
        else
        macro_changing_value = 'no_macro_chang';
        end 
        unix(['cp /home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/createGidFigures3arguments.tcl ',Path_micros,'createGidFigures3arguments.tcl'])
        unix(['/opt/GiDx64/12.1.2d/gid_offscreen -offscreen -t "source ',file_tcl_name,'"']);
        crop_images_micro(micros_evol,[Path_micros,num2str(imicro)],macro_changing_value);
        unix(['rm ',file_tcl_name]);
    end
end

function write_tcl_file(file_tcl_name,Path,file,nstep,case_post,macro_changing_value)

file_micro = [];
for iter = 1:nstep
   file_micro = [file_micro,'"',file,'_',num2str(iter),'.res" '];
end

fid = fopen(file_tcl_name,'w+');
fprintf(fid,'GiD_Process PostProcess \n');
fprintf(fid,['set Path_video_output "',file,'"\n']);
fprintf(fid,['set Path_micros_input [list ',file_micro,']\n']);

fprintf(fid,['set Name_figure ',macro_changing_value,'\n']);


switch case_post
    case 'video'    
        fprintf(fid,['source "',Path,'createGidVideo.tcl"\n']);
    case 'images'
        fprintf(fid,['source "',Path,'createGidFigures3arguments.tcl"\n']);
end
deformation = 0;
switch file
    case [Path,'Perfil_Alar']
     deformation = 0.0005;   
    otherwise
     deformation = 0.023;   
end

fprintf(fid,['VideoPrueba ',num2str(deformation) ,' $Path_video_output $Path_micros_input $Name_figure\n']);
fprintf(fid,['GiD_Process Mescape Quit']);
fclose(fid);

end

function crop_images_micro(micros_evol,Path,macro_changing_value)
for istep = 1:size(micros_evol,2)
    name_file = [Path,'-',num2str(istep),'_',macro_changing_value,'.png'];
    unix(['convert -crop 500x500+0+0 -gravity Center ',name_file,' ',name_file]);
    
end

end