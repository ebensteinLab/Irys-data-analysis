function read_result=read_write_try(type,input1,mapped_molecule,handles,settings)
err_pause=2;
path=mapped_molecule{end-1};path=path{:};
ff1=strfind(path,'\');

switch type
    case 'stitch'
        scans_directories=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                localdir=rdir([settings.directory,'FOV_files\',path,'Stitch',num2str(mapped_molecule{2},'%d'),'.fov']);
                if isempty(localdir)
                    dirall_stitch=cellfun(@(x) rdir([x,path,'**\Stitch',num2str(mapped_molecule{2},'%d'),'.fov'])',scans_directories,'UniformOutput',0);
                    dirall_stitch=[dirall_stitch{:}];
                    fov_path=dirall_stitch(1).name;
                    ff1=find(fov_path=='\',1,'last');
                    copyfile([fov_path(1:ff1),'*.fov'],[settings.directory,'FOV_files\',path],'f');
                    delete([settings.directory,'FOV_files\',path,'Stitch.fov']);
                    localdir=rdir([settings.directory,'FOV_files\',path,'Stitch',num2str(mapped_molecule{2},'%d'),'.fov']);
                end
                fileID = fopen(localdir(1).name);
                [~] = textscan(fileID,'%*[^\n]',4);
                read_result=textscan(fileID,'%f %f %f %f %f %f');
                fclose(fileID);
            catch
                set(handles.status_bar,'String',['Error loading: ',path,'Detect molecules\Stitch',num2str(mapped_molecule{2},'%d'),'.fov']);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'metadata'
        scans_directories=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                if numel(ff1)>1
                    dirall=cellfun(@(x) rdir([x,path,path(ff1(1)+1:ff1(2)-1),' Metadata.xml'])',scans_directories,'UniformOutput',0);
                else
                    dirall=cellfun(@(x) rdir([x,path,path(1:end-1),' Metadata.xml'])',scans_directories,'UniformOutput',0);
                end
                dirall=[dirall{:}];
                read_result=fileread(dirall(1).name);
            catch
                if numel(ff1)>1
                    set(handles.status_bar,'String',['Error loading: ',path,path(ff1(1)+1:ff1(2)-1),' Metadata.xml']);
                else
                    set(handles.status_bar,'String',['Error loading: ',path,path(1:end-1),' Metadata.xml']);
                end
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'dir_tiff'
        scans_directories=input1;
        dirall=[];
        while isempty(dirall)
            dirall=cellfun(@(x) rdir([x,path,'*',num2str(mapped_molecule{2},'%02d'),'.tiff'])',scans_directories,'UniformOutput',0);
            dirall=[dirall{:}];
            if isempty(dirall)
                set(handles.status_bar,'String',['Error loading: ',path,'*',num2str(mapped_molecule{2},'%02d'),'.tiff']);
                pause(err_pause);
            end
        end
        read_result=dirall;
        
    case 'creat_dir'
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                mkdir([settings.directory,'background_images\']);
            catch
                set(handles.status_bar,'String',['Error creating',settings.directory,'background_images\']);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'read_backgrounds'
        bckgrnd_path=input1{1};
        dir_bckgrnd=input1{2};
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                read_result(:,:,1) = imread(bckgrnd_path,1);
                read_result(:,:,2) = imread(bckgrnd_path,2);
                read_result(:,:,3) = imread(bckgrnd_path,3);
            catch
                set(handles.status_bar,'String',['Error loading: ',dir_bckgrnd(1).name]);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end  
                
    case 'write_backgrounds'
        if numel(ff1)>1
            bckgrnd_path=[settings.directory,'background_images\',path(ff1(1)+1:ff1(2)-1),'_bckgrnd',num2str(mapped_molecule{2},'%02d'),'.tiff'];
        else
            bckgrnd_path=[settings.directory,'background_images\',path(1:end-1),'_bckgrnd',num2str(mapped_molecule{2},'%02d'),'.tiff'];
        end
        Bbackground=input1(:,:,1);
        Gbackground=input1(:,:,2);
        Rbackground=input1(:,:,3);
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                imwrite(Bbackground,bckgrnd_path)
                imwrite(Gbackground,bckgrnd_path,'WriteMode','append')
                imwrite(Rbackground,bckgrnd_path,'WriteMode','append')
            catch
                set(handles.status_bar,'String',['Error writing: ',bckgrnd_path]);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end

    case 'read_images'
        dirall=input1{1};
        c=input1{2};
        try_count=0;
        err_count=0;       
        while try_count == err_count
            try
                read_result(:,:,:) = tiffread2(dirall(1).name,c);%imread(dirall(1).name,c(row))%c(row)+1%c(row)+2
            catch
                set(handles.status_bar,'String',['Error reading: ',dirall(1).name]);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end  
        
    case 'write_initial_images'
        molecule_image=input1{1}(:,:,:,1);
        molecule_image_original=input1{1}(:,:,:,2);
        molecule_image_background=input1{1}(:,:,:,3);
        L=input1{2};
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                writing_initial_image_files(L,molecule_image,'molecules_images',mapped_molecule,settings);
                writing_initial_image_files([],molecule_image_original,'molecules_images_originals',mapped_molecule,settings);
                writing_initial_image_files([],molecule_image_background,'molecules_images_backgrounds',mapped_molecule,settings);
            catch
                set(handles.status_bar,'String','Error writing images');
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'load_pre_images'
        dirall_pre_image_save3=input1{1};
        dirall_pre_image_save2=input1{2};
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                read_result(:,:,:,1)=imread(dirall_pre_image_save3(1).name);
                read_result(:,:,:,2)=imread(dirall_pre_image_save2(1).name);
            catch
                set(handles.status_bar,'String',['Error reading: ',dirall_pre_image_save3(1).name,' and ',dirall_pre_image_save2(1).name]);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'load_final_image'
        dirall_pre_image_save1=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                read_result=imread(dirall_pre_image_save1(1).name);
            catch
                set(handles.status_bar,'String',['Error reading: ',dirall_pre_image_save1(1).name]);
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'save_image'
        molecule_image=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                writing_initial_image_files([],molecule_image,'molecules_images',mapped_molecule,settings);
            catch
                set(handles.status_bar,'String','Error writing images');
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'write_nonmapped_bedgraphs'
        Rstrand=input1{1};
        Gstrand=input1{2};        
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                runtime=mapped_molecule{19};
                runtime=runtime{:};
                ff4=strfind(runtime,'_');
                runtime=runtime(ff4(end-2)+1:end-1);
                if ~exist([settings.directory,'bedgraphs\'], 'dir')
                    mkdir([settings.directory,'bedgraphs\']);
                end
                savefilename=[settings.directory,'bedgraphs\','run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(mapped_molecule{1},'%d'),'_length',num2str(round(mapped_molecule{10}),'%dPx'),'_ch2.txt'];
                fileID1 = fopen(savefilename,'w');                
                fprintf(fileID1,'%d\r\n',round(Gstrand));
                fclose(fileID1);
                savefilename=[settings.directory,'bedgraphs\','run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(mapped_molecule{1},'%d'),'_length',num2str(round(mapped_molecule{10}),'%dPx'),'_ch1.txt'];
                fileID2 = fopen(savefilename,'w');                
                fprintf(fileID2,'%d\r\n',round(Rstrand));
                fclose(fileID2);                
            catch
                set(handles.status_bar,'String','Error writing nonmapped bedgraphs');
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'database_mol'
        dirall_mol=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                fileID = fopen(dirall_mol);
                [~] = textscan(fileID,'%*[^\n]',4);
                read_result=textscan(fileID,'%d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %d %f %f %d %d');
                fclose(fileID);
            catch
                set(handles.status_bar,'String','Error loading mol file');
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'database_bnx'
        dirall=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                dirall2=rdir([dirall.name(1:find(dirall.name=='\',1,'last')),'*.bnx']);
                read_result=fileread(dirall2(1).name);
            catch
                set(handles.status_bar,'String','Error loading bnx file');
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
   
    case 'database_xmap'
        dirall=input1;
        try_count=0;
        err_count=0;
        while try_count == err_count
            try
                fileID = fopen(dirall.name);
                found=false;
                while ~found
                    readline=fgetl(fileID);
                    found=~isempty(strfind(readline,'#f'));
                end
                readline=readline(3:end);
                readline=strrep(readline, 'int', '%d');
                readline=strrep(readline, 'float', '%f');
                readline=strrep(readline, 'string', '%s');
                ff1=strfind(readline,'%');
                readline(ff1(2)+1)='s';
                read_result = textscan(fileID,readline);
                fclose(fileID);
            catch
                set(handles.status_bar,'String','Error loading xmap file');
                err_count=err_count+1;
                pause(err_pause);
            end
            try_count=try_count+1;
        end
        
    case 'write_mol_length_diff'
        %{
        savefilename=[settings.directory,'mol_lengths_diffs.txt'];
        if isempty(dir(savefilename))
            fileID1 = fopen(savefilename,'w');
        else
            fileID1 = fopen(savefilename,'a');
        end
        fprintf(fileID1,[num2str(abs(calc_mol_length-double(mapped_molecule{10})),'%3.2f'),'\r\n']);
        fclose(fileID1);
        %}
end