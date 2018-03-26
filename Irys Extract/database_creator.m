function database_creator(settings,handles)
if ~exist([settings.directory,'database\'], 'dir')
    mkdir([settings.directory,'database\']);
end
%for saving loading time, check if theres a presaved list of *.mol files.
%if there is one then load it. if not then create it.
if isempty(dir([settings.directory,'database\dirall_mol.mat'])) || settings.build_mol_list
    dirall_mol_org=cellfun(@(x) rdir([x,'**\Detect Molecules\molecules*.mol'])',settings.scans_directories,'UniformOutput',0);
    dirall_mol_org=[dirall_mol_org{:}];
    save([settings.directory,'database\dirall_mol.mat'],'dirall_mol_org');
else
    load([settings.directory,'database\dirall_mol.mat']);
end

save_length=5e5;
count=0;
count2=0;

bnx_flag=0;
dirall=rdir([settings.xmap_directory,'**\MoleculeQualityReport.xmap']);
if isempty(dirall)
    dirall=rdir([settings.xmap_directory,'**\Molecules.bnx']);
    bnx_flag=1;
    if isempty(dirall)
        set(handles.status_bar,'String','Can''t find xmap or bnx files');
        drawnow;
    end
end

mapped_molecules_all=cell(1,20);
mapped_molecules_all(:) = {[]};
nonmapped_molecules_all=cell(1,20);
nonmapped_molecules_all(:) = {[]};
for i2=1:length(dirall)
    mapped_molecules_MQR=cell(1,20);
    mapped_molecules_MQR(:) = {[]};
    
    if ~bnx_flag
        %read the xmap file
        set(handles.status_bar,'String',['Reading xmap ',num2str(i2,'%d'),'/',num2str(length(dirall),'%d')]);
        drawnow;
        C=read_write_try('database_xmap',dirall(i2),{{[]},{[]}},[],[]);
        
        %create a list of the mapped mol IDs
        long_molIDs=C{2};
        if length(long_molIDs{1})>18
            molIDs=cellfun(@(x) x(14:end),long_molIDs,'UniformOutput',0);
            molIDs=cellfun(@str2num,molIDs,'UniformOutput',0);
            molIDs=cell2mat(molIDs);
            molIDs=uint32(molIDs);
        else
            molIDs=cellfun(@str2num,long_molIDs,'UniformOutput',0);
            molIDs=cell2mat(molIDs);
        end
    else
        C=cell(1,15);
        molIDs=[];
    end
    
    set(handles.status_bar,'String',['Reading bnx ',num2str(i2,'%d'),'/',num2str(length(dirall),'%d')]);
    drawnow;
    file_text=read_write_try('database_bnx',dirall(i2),{{[]},{[]}},[],[]);
    set(handles.status_bar,'String',['Processing bnx ',num2str(i2,'%d'),'/',num2str(length(dirall),'%d')]);
    drawnow;
    
    if numel(file_text)>1000000
        text_length=1000000;
    else
        text_length=numel(file_text);
    end
    ff3=strfind(file_text(1:text_length), [sprintf('\n'),'# Run Data']);
    ff4=strfind(file_text(ff3(end):ff3(end)+1000),sprintf('\n'));
    ff3=[ff3,ff3(end)+ff4(2)-1];
    runIDs=cell(length(ff3)-1,1);
    for i4=1:length(ff3)-1
        readline=file_text(ff3(i4)+1:ff3(i4+1)-1);
        ff2=find(readline=='_',3,'last');
        runIDs{i4}=readline(ff2(1)+1:ff2(3)+2);
    end
    ff3=strfind(file_text(1:text_length), [sprintf('\n'),'0']);
    file_text=file_text(ff3(1):end); text_length=text_length-ff3(1)+1;
    ff3=strfind(file_text,sprintf('\n'));
    ff4=ff3(2:end);
    ff3=ff3(1:end-1);
    count1=0;
    if isempty(strfind(file_text(1:text_length),'QX22'))
        line_skip=4;
    else
        line_skip=7;
    end        
    bnx0=cell(length(ff3)./line_skip,1);
    tic
    for i4=1:line_skip:length(ff3)
        count1=count1+1;
        bnx0{count1}=file_text(ff3(i4):ff4(i4));
    end
    bnx0=cellfun(@(x) textscan(x, '%f %f %f %f %f %f %f %f %f %s %f %f %f \n', 1),bnx0,'UniformOutput',0);
    bnx0=cellfun(@(x) [x(1:9),x(11:13)],bnx0,'UniformOutput',0);
    bnx0=cell2mat(cellfun(@(x) cell2mat(x),bnx0,'UniformOutput',0));
    bnx1=zeros(max(bnx0(:,2)),12);
    bnx1(bnx0(:,2),:)=bnx0;
    bnx0=bnx1;clear bnx1;
    
    D{1}=bnx0(:,7);%OriginalMoleculeId
    D{15}=bnx0(:,8);%ScanNumber
    D{2}=bnx0(:,11);%RunId
    D{16}=(1:size(bnx0,1))';%max(molIDs);
    for i5=3:14
        D{i5}(molIDs)=C{i5};
        if i5==8 || i5==10 || i5==14
            D{i5}(size(bnx0,1)-1:size(bnx0,1))=cell(2,1);
        else
            D{i5}(size(bnx0,1)-1:size(bnx0,1))=[0 0];
        end
        D{i5}=D{i5}';
    end
    C=D;
    clear D;
        
    [~,sort_ind_scans] = sort(C{15},'ascend');
    sort_ind_scans_cell=cell(1,numel(C));
    sort_ind_scans_cell(:) = {sort_ind_scans};
    C=cellfun(@(x,y) x(y),C,sort_ind_scans_cell,'UniformOutput',0);
    
    set(handles.status_bar,'String',['Done processing bnx ',num2str(i2,'%d'),'/',num2str(length(dirall),'%d')]);
    drawnow;
    
    for i1=1:length(runIDs)
        %only xmap entries from the current run
        current_runID_cell=cell(1,length(C));
        current_runID_cell(:) = {find(C{2}==i1)};
        current_runID_C=cellfun(@(x1,x2) x1(x2),C,current_runID_cell,'UniformOutput',0);
        
        for i3=1:max(current_runID_C{15})
            set(handles.status_bar,'String',['bnx ',num2str(i2,'%d'),'(',num2str(length(dirall),'%d'),...
                '), run ',num2str(i1,'%d'),'(',num2str(length(runIDs),'%d'),...
                '), scan ',num2str(i3,'%d'),'(',num2str(max(current_runID_C{15}),'%d'),')']);
            drawnow;
            %search for the correct mol file path
            dirall_mol={dirall_mol_org.name};
            dataset_cell=cell(1,length(dirall_mol));
            dataset_cell(:) = {runIDs{i1}};
            dirall_mol=dirall_mol(cellfun(@isempty,cellfun(@strfind,dirall_mol,dataset_cell,'UniformOutput',0))==0);
            dataset_cell=cell(1,length(dirall_mol));            
            dataset_cell(:) = {['Molecules',num2str(i3),'.mol']};
            dirall_mol=dirall_mol(cellfun(@isempty,cellfun(@strfind,dirall_mol,dataset_cell,'UniformOutput',0))==0);
            if length(dirall_mol)~=1
                set(handles.status_bar,'String',['Error loading mol file from run ',runIDs{i1}]);
                drawnow;
                break;
            end
            dirall_mol=dirall_mol{1};

            %load the mol file
            C2=read_write_try('database_mol',dirall_mol,{{[]},{[]}},[],[]);
            
            %only xmap entries from the current scan
            current_runID_cell=cell(1,length(current_runID_C));
            current_runID_cell(:) = {find(current_runID_C{15}==i3)};
            current_runID_C2=cellfun(@(x1,x2) x1(x2),current_runID_C,current_runID_cell,'UniformOutput',0);
            
            current_molID_cell=cell(1,length(C2));
            current_molID_cell(:) = {current_runID_C2{1}};
            relevant_C2_cell=cellfun(@(x1,x2) x1(x2),C2,current_molID_cell,'UniformOutput',0);
            relevant_C2_cell([1,2,3,5,13:16,18:23])=[];
            %first one is molIDs from the xmap, after that are the remaining mol entries
            %and after that the data from the xmap with some modification
            %from the bnx
            mapped_molecules=[current_runID_C2(16),relevant_C2_cell,current_runID_C2([3:9,14])];
            
            %add a column of the path of the dataset
            data_path=dirall_mol;
            for i7=1:length(settings.scans_directories)
                temp_scan_dir=settings.scans_directories{i7};
                if numel(data_path)>numel(temp_scan_dir)
                    if strcmp(data_path(1:length(temp_scan_dir)),temp_scan_dir)
                        ff4=find(data_path=='\',2,'last');
                        data_path=data_path(length(temp_scan_dir)+1:ff4(1));
                        break;
                    end
                else
                    continue;
                end
            end
            %ff4=find(data_path=='\',4,'last');
            %data_path=data_path(ff4(1)+1:ff4(3));
            data_path_cell=cell(length(mapped_molecules{1}),1);
            data_path_cell(:) = {data_path};
            mapped_molecules{numel(mapped_molecules)+1}=data_path_cell;
            
            %add a column of the path of the xmap
            data_path=dirall(i2).name;
            ff4=find(data_path=='\',2,'last');
            %data_path=data_path(length(settings.xmap_directory)+1:ff4);
            data_path=data_path(ff4(1)+1:ff4(2)-1);
            data_path_cell=cell(length(mapped_molecules{1}),1);
            data_path_cell(:) = {data_path};
            mapped_molecules{numel(mapped_molecules)+1}=data_path_cell;
            
            mapped_molecules_MQR=cellfun(@(x1,x2) [x1;x2],mapped_molecules_MQR,mapped_molecules,'UniformOutput',0);
        end
    end
    
    ff1=mapped_molecules_MQR{17}>0;
    nonmapped_molecules_MQR=cellfun(@(x) x(~ff1),mapped_molecules_MQR,'UniformOutput',0);
    mapped_molecules_MQR=cellfun(@(x) x(ff1),mapped_molecules_MQR,'UniformOutput',0);
    
    confidences=mapped_molecules_MQR{17};
    [~,sort_ind_conf] = sort(confidences,'descend');
    sort_ind_conf_cell=cell(1,numel(mapped_molecules_MQR));
    sort_ind_conf_cell(:) = {sort_ind_conf};
    mapped_molecules_MQR=cellfun(@(x,y) x(y),mapped_molecules_MQR,sort_ind_conf_cell,'UniformOutput',0);
    
    mapped_molecules_all=cellfun(@(x1,x2) [x1;x2],mapped_molecules_all,mapped_molecules_MQR,'UniformOutput',0);
    nonmapped_molecules_all=cellfun(@(x1,x2) [x1;x2],nonmapped_molecules_all,nonmapped_molecules_MQR,'UniformOutput',0);
    while length(mapped_molecules_all{1})>=save_length
        count=count+1;
        mapped_molecules_all_temp=cellfun(@(x) x(1:save_length),mapped_molecules_all,'UniformOutput',0);
        save([settings.directory,'database\mapped_molecules',num2str(count,'%d'),'.mat'],'mapped_molecules_all_temp');
        if length(mapped_molecules_all{1})>save_length
            mapped_molecules_all=cellfun(@(x) x(save_length+1:end),mapped_molecules_all,'UniformOutput',0);
        else
            mapped_molecules_all(:) = {[]};
        end
    end
    while length(nonmapped_molecules_all{1})>=save_length
        count2=count2+1;
        mapped_molecules_all_temp=cellfun(@(x) x(1:save_length),nonmapped_molecules_all,'UniformOutput',0);
        save([settings.directory,'database\nonmapped_molecules',num2str(count2,'%d'),'.mat'],'mapped_molecules_all_temp');
        if length(nonmapped_molecules_all{1})>save_length
            nonmapped_molecules_all=cellfun(@(x) x(save_length+1:end),nonmapped_molecules_all,'UniformOutput',0);
        else
            nonmapped_molecules_all(:) = {[]};
        end
    end  
end
if ~isempty(mapped_molecules_all{1})
    count=count+1;
    mapped_molecules_all_temp=mapped_molecules_all;
    save([settings.directory,'database\mapped_molecules',num2str(count,'%d'),'.mat'],'mapped_molecules_all_temp');
end
if ~isempty(nonmapped_molecules_all{1})
    count2=count2+1;
    mapped_molecules_all_temp=nonmapped_molecules_all;
    save([settings.directory,'database\nonmapped_molecules',num2str(count2,'%d'),'.mat'],'mapped_molecules_all_temp');
end