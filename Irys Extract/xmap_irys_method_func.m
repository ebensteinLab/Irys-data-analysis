function xmap_irys_method_func(settings,handles)
global i3 mapping_limits mapping_summary averaged_molecules;
directory=settings.directory;
interest_zones=settings.interest_zones;
bedres=settings.bed_res;
chr_lengths=settings.chr_lengths;

if isempty(dir([directory,'database\nonmapped_molecules*.mat']))
    database_creator(settings,handles);
end
mapped_molecules_all=cell(1,20);
mapped_molecules_all(:) = {[]};
dirall=rdir([directory,'database\mapped_molecules*.mat']);
for i1=1:length(dirall)
    load(dirall(i1).name);
    mapped_molecules_all=cellfun(@(x1,x2) [x1;x2],mapped_molecules_all,mapped_molecules_all_temp,'UniformOutput',0);
end
if settings.include_non_mapped
    dirall=rdir([directory,'database\nonmapped_molecules*.mat']);
    for i1=1:length(dirall)
        load(dirall(i1).name);
        mapped_molecules_all=cellfun(@(x1,x2) [x1;x2],mapped_molecules_all,mapped_molecules_all_temp,'UniformOutput',0);
    end
end
clear mapped_molecules_all_temp;

path_old=[];
total_molecules_number=length(mapped_molecules_all{1});

if settings.include_non_mapped
    ff1=mapped_molecules_all{10}*440>=settings.min_mol_length*1000;
    mapped_molecules_all=cellfun(@(x) x(ff1),mapped_molecules_all,'UniformOutput',0);
else
    ff1=mapped_molecules_all{17}>=settings.confidence_lower_limit;
    mapped_molecules_all=cellfun(@(x) x(ff1),mapped_molecules_all,'UniformOutput',0);
    
    ff1=abs(double(mapped_molecules_all{15})-double(mapped_molecules_all{14}))>=settings.min_mol_length*1000;
    mapped_molecules_all=cellfun(@(x) x(ff1),mapped_molecules_all,'UniformOutput',0);
    
    if settings.interest_zones_check
        ff1=logical(zeros(numel(mapped_molecules_all{17}),1));
        for i11=1:length(interest_zones)
            interest_zone_attr=sscanf(interest_zones{i11}, '%d,%d,%d');
            interest_chr=interest_zone_attr(1);
            interest_zone_top=interest_zone_attr(2);
            interest_zone_bottom=interest_zone_attr(3);
            ff1=ff1 | (mapped_molecules_all{11}==interest_chr & ((mapped_molecules_all{14}<=interest_zone_top & mapped_molecules_all{15}>=interest_zone_top) | (mapped_molecules_all{14}<=interest_zone_bottom & mapped_molecules_all{15}>=interest_zone_bottom)  | (mapped_molecules_all{14}>=interest_zone_bottom & mapped_molecules_all{15}<=interest_zone_top)  | (mapped_molecules_all{14}<=interest_zone_bottom & mapped_molecules_all{15}>=interest_zone_top)));
        end
        mapped_molecules_all=cellfun(@(x) x(ff1),mapped_molecules_all,'UniformOutput',0);
    end
end

if settings.selected_IDs_check
    selected_IDs_list_inds=[0,strfind(settings.selected_IDs_list,','),length(settings.selected_IDs_list)+1];
    ff1=[];
    for i5=1:length(selected_IDs_list_inds)-1
        ff1=[ff1,find(mapped_molecules_all{1}==str2num(settings.selected_IDs_list(selected_IDs_list_inds(i5)+1:selected_IDs_list_inds(i5+1)-1)))];
    end
    mapped_molecules_all=cellfun(@(x) x(ff1),mapped_molecules_all,'UniformOutput',0);
end

%ff1=cellfun(@(x) ~isempty(x),strfind(mapped_molecules_all{19},'2016-04'));
%mapped_molecules_all=cellfun(@(x) x(ff1),mapped_molecules_all,'UniformOutput',0);

filtered_molecules_number=length(mapped_molecules_all{1});

if ~exist([directory,'autosave\'], 'dir')
    mkdir([directory,'autosave\']);
end

if ~settings.continue_past
    delete([directory,'bed_channel*.txt']);
    delete([directory,'bedgraphs_reduced\*.txt']);
    delete([directory,'autosave\*.txt']);
    mapping_limits=cell(1,0);
    if settings.write_bed
        mapping_summary=cell(1,24);
    end
    averaged_molecules=cell(1,24);
    for i=1:24
        averaged_molecules{i}=zeros(ceil(chr_lengths(i)./bedres),3);
    end
    past_run_i3=1;
    set(handles.hmc_marked_stat,'String','0');
    set(handles.mapping_OK_stat,'String','0');
    set(handles.range_OK_stat,'String','0');
    set(handles.failed_fit_stat,'String','0');
    set(handles.no_peaks_stat,'String','0');
else
    load([directory,'past_run_data.mat']);
    set(handles.hmc_marked_stat,'String',num2str(status_counters(7),'%d'));
    set(handles.mapping_OK_stat,'String',num2str(status_counters(3),'%d'));
    set(handles.range_OK_stat,'String',num2str(status_counters(4),'%d'));
    set(handles.failed_fit_stat,'String',num2str(status_counters(5),'%d'));
    set(handles.no_peaks_stat,'String',num2str(status_counters(6),'%d'));
end
dir_molecules=dir([settings.directory,'molecules_images\*.tiff']);
dir_molecules2=cell(length(dir_molecules),1);
for i1=1:length(dir_molecules)
    dir_molecules2{i1}=dir_molecules(i1).name;
end
dir_molecules=strjoin(dir_molecules2');
clear dir_molecules2;
tic
        
for i3=past_run_i3:filtered_molecules_number
    if logical(get(handles.save_stop,'Value'))
        past_run_i3=i3;
        status_counters(7)=str2double(get(handles.hmc_marked_stat,'String'));
        status_counters(3)=str2double(get(handles.mapping_OK_stat,'String'));
        status_counters(4)=str2double(get(handles.range_OK_stat,'String'));
        status_counters(5)=str2double(get(handles.failed_fit_stat,'String'));
        status_counters(6)=str2double(get(handles.no_peaks_stat,'String'));
        save([directory,'past_run_data.mat'],'past_run_i3','mapping_limits','mapping_summary','averaged_molecules','status_counters');
        i3=i3-1;
        break;
    end
    if toc>settings.autosave_time*60
        past_run_i3=i3;
        status_counters(7)=str2double(get(handles.hmc_marked_stat,'String'));
        status_counters(3)=str2double(get(handles.mapping_OK_stat,'String'));
        status_counters(4)=str2double(get(handles.range_OK_stat,'String'));
        status_counters(5)=str2double(get(handles.failed_fit_stat,'String'));
        status_counters(6)=str2double(get(handles.no_peaks_stat,'String'));
        save([directory,'past_run_data.mat'],'past_run_i3','mapping_limits','mapping_summary','averaged_molecules','status_counters');
        if ~isempty(dir([settings.directory,'bedgraphs_reduced\*.txt']))
            copyfile([settings.directory,'bedgraphs_reduced\*.txt'],[settings.directory,'autosave\'],'f');
        end
        tic
    end
    set(handles.status_bar,'String',[num2str(i3,'%d'),' / ',num2str(filtered_molecules_number,'%d'),' (',num2str(total_molecules_number,'%d'),' total)...']);
    mapped_molecule=cellfun(@(x) x(i3),mapped_molecules_all,'UniformOutput',0);
    mapped_molecule=mapped_molecule(:);
    set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'Producing image...']); 
    drawnow;
    
    %{
    %for extracting a list of the best interest zones
    devided_x_ind=ceil(min([double(mapped_molecule{14}),double(mapped_molecule{15})])/bedres):floor(max([double(mapped_molecule{14}),double(mapped_molecule{15})])/bedres);
    if length(averaged_molecules{mapped_molecule{11}}(:,3))>=max(devided_x_ind)
        averaged_molecules{mapped_molecule{11}}(devided_x_ind,3)=averaged_molecules{mapped_molecule{11}}(devided_x_ind,3)+ones(length(devided_x_ind),1);
    end
    continue;
    %}
    
    MQR=mapped_molecule{20};
    MQR=MQR{:};    
    if ~isempty(strfind(dir_molecules,['_ID',num2str(mapped_molecule{1},'%d'),'_MQR',MQR])) && ~(settings.write_bed || settings.write_graph)
        continue;
    end
    
    [~,Rstrand,Gstrand]=producing_molecule_image_and_intensity_profiles(mapped_molecule,settings,handles);
    
    if logical(get(handles.save_stop,'Value'))
        past_run_i3=i3;
        status_counters(7)=str2double(get(handles.hmc_marked_stat,'String'));
        status_counters(3)=str2double(get(handles.mapping_OK_stat,'String'));
        status_counters(4)=str2double(get(handles.range_OK_stat,'String'));
        status_counters(5)=str2double(get(handles.failed_fit_stat,'String'));
        status_counters(6)=str2double(get(handles.no_peaks_stat,'String'));
        save([directory,'past_run_data.mat'],'past_run_i3','mapping_limits','mapping_summary','averaged_molecules','status_counters');
        break;
    end
    
    [maxtab,~]=peakdet(Gstrand,settings.hmc_threshold,1:length(Gstrand));
    if numel(maxtab)/2<settings.hmc_num_peaks
        continue;
    end
    set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'2nd channel OK...']);
    set(handles.hmc_marked_stat,'String',num2str(str2double(get(handles.hmc_marked_stat,'String'))+1,'%d'));
    
    if ~(settings.write_bed || settings.write_graph)
        continue;
    end
    
    align_temp=mapped_molecule{18};    
    if isempty(align_temp{:})
        read_write_try('write_nonmapped_bedgraphs',{Rstrand,Gstrand},mapped_molecule,handles,settings);
        continue;
    end
    
    set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'Fitting molecule...']);
    drawnow;
    
    if ~exist('qcmap', 'var')
        qcmap=[];rcmap=[];
    end
    [modified_x,Rstrand,Gstrand,path_old,qcmap,rcmap]=genome_alignment(path_old,mapped_molecule,Rstrand,Gstrand,qcmap,rcmap,settings,handles);
    
    if isempty(modified_x)
        continue;
    end
    set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'alignment OK...']);
    set(handles.mapping_OK_stat,'String',num2str(str2double(get(handles.mapping_OK_stat,'String'))+1,'%d'));
      
    if mapped_molecule{11}>24
        mapped_molecule{11}=1;
        modified_x=modified_x+4e6;
    end
    if modified_x(1)<0 || modified_x(end)>chr_lengths(mapped_molecule{11}) || (modified_x(end)-modified_x(1))/(double(mapped_molecule{15})-double(mapped_molecule{14}))>1+settings.mapping_range_extra/100
        continue;
    end
    
    set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'mapping range OK']);
    set(handles.range_OK_stat,'String',num2str(str2double(get(handles.range_OK_stat,'String'))+1,'%d'));
    set(get(handles.axes1,'Title'),'String',[get(get(handles.axes1,'Title'),'String'),', mapping range OK']);
    drawnow;
    
    if settings.write_bed
        mapping_summary{mapped_molecule{11}}=[mapping_summary{mapped_molecule{11}};[modified_x(1),modified_x(end)]];%[mapped_molecule{14},mapped_molecule{15}]];
    end
    
    mapping_limits=writing_bed_and_image_files(mapping_limits,mapped_molecule,modified_x,Gstrand,Rstrand,settings);
    
    devided_x_ind=ceil(modified_x./bedres);
    averaged_molecules{mapped_molecule{11}}(devided_x_ind,1)=averaged_molecules{mapped_molecule{11}}(devided_x_ind,1)+Rstrand;
    averaged_molecules{mapped_molecule{11}}(devided_x_ind,2)=averaged_molecules{mapped_molecule{11}}(devided_x_ind,2)+Gstrand;
    averaged_molecules{mapped_molecule{11}}(devided_x_ind,3)=averaged_molecules{mapped_molecule{11}}(devided_x_ind,3)+ones(length(devided_x_ind),1);
end

%{
set(handles.status_bar,'String','Writing region file...');
drawnow;
if write_bed
    write_bed_region_file(directory,mapping_summary);
end
%}
set(handles.status_bar,'String','Calculating and writing consensus file...');
drawnow;
calculate_and_write_consensus(directory,averaged_molecules,bedres,chr_lengths);
if isempty(i3)
    i3=0;
end
set(handles.status_bar,'String',['Done processing ',num2str(i3,'%d'),' molecules out of ',num2str(filtered_molecules_number,'%d'),' from a total of ',num2str(total_molecules_number,'%d')]);