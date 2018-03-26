function writing_initial_image_files(L,molecule_image,folder_name,mapped_molecule,settings);

if ~exist([settings.directory,folder_name,'\'], 'dir')
    mkdir([settings.directory,folder_name,'\']);
end
runtime=mapped_molecule{19};
runtime=runtime{:};
ff4=strfind(runtime,'_');
runtime=runtime(ff4(end-2)+1:end-1);
MQR=mapped_molecule{20};
MQR=MQR{:};
%count=1;
%savefilename=[settings.directory,folder_name,'\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(mapped_molecule{1},'%d'),'_length',num2str(round((abs(double(mapped_molecule{15})-double(mapped_molecule{14})))/1000),'%dkbp'),'_not_mapped',num2str(count,'%d'),'mol.tiff'];
%while exist(savefilename, 'file')
%    count=count+1;
%    savefilename=[settings.directory,folder_name,'\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(mapped_molecule{1},'%d'),'_length',num2str(round((abs(double(mapped_molecule{15})-double(mapped_molecule{14})))/1000),'%dkbp'),'_not_mapped',num2str(count,'%d'),'mol.tiff'];
%end

align_temp=mapped_molecule{18};
if isempty(align_temp{:})
    savefilename=[settings.directory,folder_name,'\','run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(mapped_molecule{1},'%d'),'_MQR',MQR,'_length',num2str(round(mapped_molecule{10}),'%dPx'),'.tiff'];       
else
    savefilename=[settings.directory,folder_name,'\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(mapped_molecule{1},'%d'),'_MQR',MQR,'_length',num2str(round((abs(double(mapped_molecule{15})-double(mapped_molecule{14})))/1000),'%dkbp'),'_not_mapped.tiff'];
end
imwrite(molecule_image,savefilename);

if settings.write_images_FOV && ~isempty(L)
    savefilename=[savefilename(1:end-5),'FOVs.tiff'];
    imwrite(L,savefilename);
end