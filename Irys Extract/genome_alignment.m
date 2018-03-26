function [modified_x,Rstrand,Gstrand,path_old,qcmap,rcmap]=genome_alignment(path_old,mapped_molecule,Rstrand,Gstrand,qcmap,rcmap,settings,handles)

path=mapped_molecule{end};path=path{:};
%path=[settings.xmap_directory,path];
if ~strcmp(path,path_old)
    path1=[settings.xmap_directory,path,'\'];
    dirall=rdir([path1,'*_q.cmap']);
    if isempty(dirall)
        path1=settings.xmap_directory;
        dirall=rdir([path1,'*_q.cmap']);
    end
    fileID = fopen(dirall(1).name);    
    found=false;
    while ~found
        readline=fgetl(fileID);
        found=~isempty(strfind(readline,'#f'));
    end
    readline=readline(3:end);
    readline=strrep(readline, 'int', '%f');
    readline=strrep(readline, 'float', '%f');    
    %[~] = textscan(fileID,'%*[^\n]',11);%was 9 lines % was 10 lines
    qcmap=textscan(fileID,readline);%was 9 values
    fclose(fileID);
    qcmap=cell2mat(qcmap);
        
    dirall=rdir([path1,'*_r.cmap']);    
    fileID = fopen(dirall(1).name);
    found=false;
    while ~found
        readline=fgetl(fileID);
        found=~isempty(strfind(readline,'#f'));
    end
    readline=readline(3:end);
    readline=strrep(readline, 'int', '%f');
    readline=strrep(readline, 'float', '%f');    
    %[~] = textscan(fileID,'%*[^\n]',11); %was 10 lines
    rcmap=textscan(fileID,readline);
    fclose(fileID);
    rcmap=cell2mat(rcmap);
    path_old=path;
end
current_qcmap=qcmap(qcmap(:,1)==double(mapped_molecule{1}),:);%mol ID filter
bspq_labels=current_qcmap(current_qcmap(:,5)==1,6); %label channel filter
bspq_labels=bspq_labels./(current_qcmap(1,2)./double(mapped_molecule{10}))+1;
alignment=mapped_molecule{18};alignment=alignment{:};
ff1=find(alignment=='(')+1;
ff2=find(alignment==',');
ff3=find(alignment==')')-1;
alignment_arr=zeros(length(ff1),2);
for i4=1:length(ff1)
    alignment_arr(i4,2)=str2double(alignment(ff1(i4):ff2(i4)-1));
    alignment_arr(i4,1)=str2double(alignment(ff2(i4)+1:ff3(i4)));
end
current_rcmap=rcmap(rcmap(:,1)==double(mapped_molecule{11}),:); %chromosome filter
alignment_arr(:,2)=current_rcmap(alignment_arr(:,2),6);%from ref tag id to base position
axes(handles.axes1);
plot(alignment_arr(:,2)./1e6,-1.*ones(length(alignment_arr(:,2)),1),'.k')
hold on;
alignment_arr(:,1)=bspq_labels(alignment_arr(:,1));%from tag id to pixel
ff4=find(diff(alignment_arr(:,1))==0);
alignment_arr(ff4,2)=(alignment_arr(ff4,2)+alignment_arr(ff4+1,2))./2;
alignment_arr(ff4+1,:)=[];
original_x=1:length(Rstrand);
map_direction=mapped_molecule{16};map_direction=map_direction{:};
if map_direction=='-'
    Rstrand=Rstrand(end:-1:1);
    Gstrand=Gstrand(end:-1:1);
    original_x=original_x(end:-1:1);
end
if settings.align_ch_num
    temp_strand=Gstrand;
    Gstrand=Rstrand;
    Rstrand=temp_strand;
end    

modified_x=interp1(alignment_arr(:,1),alignment_arr(:,2),original_x,'linear','extrap');
[maxtab,~]=peakdet(Rstrand,settings.peak_detection_threshold,original_x);

%if ~isempty(maxtab)
%    diffs=maxtab(:,1)*ones(1,size(alignment_arr,1))-ones(size(maxtab,1),1)*alignment_arr(:,2)';
%    maxtab(sum(diffs<settings.peak_max_diff,2)==0,:)=[];
%end
settings.shift_amp=8; settings.shift_res=0.2;settings.stretch_amp=5;settings.stretch_res=0.1;
if size(maxtab,1)>=settings.min_num_of_peaks && size(maxtab,1)>0
    fit_results=zeros( floor(settings.shift_amp./settings.shift_res).*2+1  ,3);
    count=0;
    for stretch=1-settings.stretch_amp/100:settings.stretch_res/100:1+settings.stretch_amp/100
        for shift=-settings.shift_amp:settings.shift_res:settings.shift_amp
            new_maxtab=shift+stretch.*maxtab(:,1);
            new_maxtab=interp1(alignment_arr(:,1),alignment_arr(:,2),new_maxtab,'linear','extrap');
            diffs=new_maxtab(:,1)*ones(1,size(alignment_arr,1))-ones(size(new_maxtab,1),1)*alignment_arr(:,2)';
            [vals,~]=min(abs(diffs'));
            vals(vals>settings.shift_amp*500)=[];
            vals(vals>mean(vals)+std(vals))=[];
            vals(vals>mean(vals)+std(vals))=[];
            count=count+1;
            fit_results(count,:)=[shift,mean(vals),stretch];
        end
    end
    [~,min_fit_ind]=min(fit_results(:,2));
    modified_x1=interp1(alignment_arr(:,1),alignment_arr(:,2),original_x.*fit_results(min_fit_ind,3)+fit_results(min_fit_ind,1),'linear','extrap');
    final_peaks_locations=interp1(alignment_arr(:,1),alignment_arr(:,2),maxtab(:,1).*fit_results(min_fit_ind,3)+fit_results(min_fit_ind,1),'linear','extrap');
    plot(modified_x./1e6,Rstrand,'b-',modified_x1./1e6,Rstrand,'g-',final_peaks_locations./1e6,maxtab(:,2),'m.')
    hold off;
    axis tight;
    if fit_results(min_fit_ind,2)<settings.good_fit_cond
        if ~isempty(dir([settings.directory,'fit_pars.mat']))
            load([settings.directory,'fit_pars.mat']);
        else
            fit_pars=[];
        end
        fit_pars=[fit_pars;fit_results(min_fit_ind,1)];
        %save([settings.directory,'fit_pars.mat'],'fit_pars');
        modified_x=modified_x1;
        %title(['Shift: ',num2str(fit_results(min_fit_ind,1),'%.1f')]);
        title('Fit OK');
    else
        modified_x=[];
        set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'failed fit']);
        set(handles.failed_fit_stat,'String',num2str(str2double(get(handles.failed_fit_stat,'String'))+1,'%d'));
    end
else
    hold off;
    modified_x=[];
    set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'not enough detectable spots']);
    set(handles.no_peaks_stat,'String',num2str(str2double(get(handles.no_peaks_stat,'String'))+1,'%d'));
end
drawnow;