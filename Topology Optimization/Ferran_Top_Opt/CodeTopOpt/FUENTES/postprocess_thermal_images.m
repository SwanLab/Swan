function postprocess_thermal_images()


typemicro = 'lamin';%'opt'; %lamin

switch typemicro
    
    case 'opt'
        Path = '/home/aferrer/Documents/Doctorat/Tesi/ElasticThermal_Regularized/VademecumThermal_optimo/';
    case 'lamin'
        Path = '/home/aferrer/Documents/Doctorat/Tesi/ElasticThermal_Regularized/VademecumThermal_LaminadoNoPeriodico/';
end

[~,num_snap] = unix(['ls ',Path,' -lR | grep ^d | wc -l']);
num_snap = str2double(num_snap) -4;

Image_Case = 'Images';
%Image_Case = 'Images4repeated';
order_index = 1:num_snap;
postprocess_thermal_dif_case(Path,num_snap,Image_Case,order_index,typemicro)


Image_Case = 'ImagesPhiOrder';
%Image_Case = 'ImagesPhiOrder4repeated';
nlevel = 7;pint = 0; [~,phi_all] = vademecum_spehere_mesh(nlevel,pint,'THERMAL');
[~,order_index] = sort(phi_all);
postprocess_thermal_dif_case(Path,num_snap,Image_Case,order_index,typemicro)




end

function  postprocess_thermal_dif_case(Path,num_snap,Image_Case,order_index,typemicro)
for  imicro = 1:num_snap
    imicro_vad = order_index(imicro);
    Path_micros = [Path,num2str(imicro_vad),'/'];
    file_tcl_name = [Path_micros,'tcl_Video_file_micro_',num2str(imicro_vad),'.tcl'];
    switch typemicro
        
        case 'opt'
            num_iter =  compute_num_lines([Path_micros,'Ch.txt']) -2;
        case 'lamin'
            num_iter = 0;
    end
    file_res = [Path_micros,'RVE04N3_',num2str(num_iter),'.flavia.res'];
    write_tcl_file(file_tcl_name,Path,file_res,imicro,Image_Case)
    %       unix(['/opt/GiDx64/12.1.2d/gid_offscreen -t "source ',file_tcl_name,'"'])
    %unix(['cp ',Path_micros,'/createGidFiguresTemplate.tcl ',Path_micros,'/createGidFigures.tcl '])
    unix(['/opt/GiDx64/12.1.2d/gid_offscreen -offscreen -t "source ',file_tcl_name,'"']);
    crop_images_micro(Path,imicro,Image_Case);
    unix(['rm ',file_tcl_name]);
end


end



function write_tcl_file(file_tcl_name,Path,file_micro,imicro,Image_Case)
fid = fopen(file_tcl_name,'w+');
fprintf(fid,'GiD_Process PostProcess \n');
fprintf(fid,['set Path_video_output "',Path,'/',Image_Case,'/',num2str(imicro),'"\n']);
fprintf(fid,['set Path_micros_input [list ',file_micro,']\n']);
%fprintf(fid,['source "',Path,'createGidFigures.tcl"\n']);
fprintf(fid,['source "',Path,'createGidFigures.tcl"\n']);
deformation = 0.023;
fprintf(fid,['VideoPrueba ',num2str(deformation) ,' $Path_video_output $Path_micros_input \n']);
% fprintf(fid,['GiD_Process Mescape Quit']);
fclose(fid);

end


function crop_images_micro(Path,imicro,Image_Case)
name_file = [Path,'/',Image_Case,'/',num2str(imicro),'.png'];
unix(['convert -crop 500x500+0+0 -gravity Center ',name_file,' ',name_file]);
end