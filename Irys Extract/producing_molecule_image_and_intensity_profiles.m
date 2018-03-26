function [molecule_image,Rstrand,Gstrand]=producing_molecule_image_and_intensity_profiles(mapped_molecule,settings,handles)
scans_directories=settings.scans_directories;
extra_length=settings.extra_length;
xedges=settings.xedges;
sum_edges=settings.sum_edges;
image_size_x=512;
image_size_y=512;
path=mapped_molecule{end-1};path=path{:};
col=double(mapped_molecule{3});

stitch=read_write_try('stitch',scans_directories,mapped_molecule,handles,settings);
stitch=cell2mat(stitch);
stitch(stitch(:,3)~=col,:)=[];

if double(mapped_molecule{7})==0 && double(mapped_molecule{4})>1
    mapped_molecule{7}=image_size_y-1+stitch(double(mapped_molecule{4}),6)-24;
    mapped_molecule{4}=double(mapped_molecule{4})-1;
end

shifts=zeros(double(mapped_molecule{5}),2);
shifts(double(mapped_molecule{4})+1:double(mapped_molecule{5}),1)=cumsum(stitch(double(mapped_molecule{4})+1:double(mapped_molecule{5}),5));%+tilt_x_shift;
shifts(double(mapped_molecule{4})+1:double(mapped_molecule{5}),2)=stitch(double(mapped_molecule{4})+1:double(mapped_molecule{5}),6);%-4;%+tilt_y_shift-4;
tilt_angle=mean(stitch(double(mapped_molecule{4}):double(mapped_molecule{5}),4));
mapped_molecule{8}=double(mapped_molecule{8})+sum(stitch(double(mapped_molecule{4})+1:double(mapped_molecule{5}),5));
sum_edges=max([sum_edges,ceil((double(mapped_molecule{8})-double(mapped_molecule{6}))/2)]);

MQR=mapped_molecule{20};
MQR=MQR{:};
savefilename=[settings.directory,'molecules_images','\','*_ID',num2str(mapped_molecule{1},'%d'),'_MQR',MQR,'*.tiff'];
dirall_pre_image_save1=rdir(savefilename);
savefilename=[settings.directory,'molecules_images_originals','\','*_ID',num2str(mapped_molecule{1},'%d'),'_MQR',MQR,'*.tiff'];
dirall_pre_image_save2=rdir(savefilename);
savefilename=[settings.directory,'molecules_images_backgrounds','\','*_ID',num2str(mapped_molecule{1},'%d'),'_MQR',MQR,'*.tiff'];
dirall_pre_image_save3=rdir(savefilename);

if isempty(dirall_pre_image_save2) || isempty(dirall_pre_image_save3) %isempty(dirall_pre_image_save1) ||
    %set(handles.status_bar,'String',[get(handles.status_bar, 'String'),'Costructing image...']);
    %drawnow;
    
    file_text=read_write_try('metadata',scans_directories,mapped_molecule,handles,settings);
    ff3=strfind(file_text, '<LasersCount>');
    ff4=strfind(file_text,'</LasersCount>');
    channels=str2double(file_text(ff3(1)+13:ff4(1)-1));%3;
    ff3=strfind(file_text, '<RelativeMagnification>');
    ff4=strfind(file_text,'</RelativeMagnification>');
    scales(3)=str2double(file_text(ff3(1)+23:ff4(1)-1));%B
    if channels==3
    scales(2)=str2double(file_text(ff3(2)+23:ff4(2)-1));%G
    scales(1)=str2double(file_text(ff3(3)+23:ff4(3)-1));%R
    else if channels==2
            scales(2)=1;%G
            scales(1)=str2double(file_text(ff3(2)+23:ff4(2)-1));%R
        end
    end
    ff3=strfind(file_text, '<ScanRowCount>');
    ff4=strfind(file_text,'</ScanRowCount>');
    rows=str2double(file_text(ff3(1)+14:ff4(1)-1));%12;
    zigzag=1;
    start_from_bottom=1;
    count_direction_arr=[1 1 -1 -1 -1 1 1 -1];
    count_direction=count_direction_arr(bin2dec(num2str([zigzag,start_from_bottom,rem(col,2)]))+1);
    if count_direction==1
        c_add=1:channels:(rows)*channels;
    else
        c_add=(rows-1)*channels+1:-channels:1;
    end
    c=c_add+(col-1).*channels.*rows;
    c=c(end:-1:1);
    
    dirall=read_write_try('dir_tiff',scans_directories,mapped_molecule,handles,settings);
    
    if ~exist([settings.directory,'background_images\'], 'dir')
        read_write_try('creat_dir',[],mapped_molecule,handles,settings);
    end
    
    ff1=strfind(path,'\');
    
    if numel(ff1)>1
        bckgrnd_path=[settings.directory,'background_images\',path(ff1(1)+1:ff1(2)-1),'_bckgrnd',num2str(mapped_molecule{2},'%02d'),'.tiff'];
    else
        bckgrnd_path=[settings.directory,'background_images\',path(1:end-1),'_bckgrnd',num2str(mapped_molecule{2},'%02d'),'.tiff'];
    end    
    dir_bckgrnd=dir(bckgrnd_path);
    if ~isempty(dir_bckgrnd)
        read_result=read_write_try('read_backgrounds',{bckgrnd_path,dir_bckgrnd},mapped_molecule,handles,settings);      
        Bbackground = read_result(:,:,1);
        Gbackground = read_result(:,:,2);
        Rbackground = read_result(:,:,3);        
    else
        background_n=settings.background_samples_per_pixel;
        backgrounds=uint32(zeros(image_size_x,image_size_y,3));
        count_background=uint16(zeros(image_size_x,image_size_y,3));
        read_result=read_write_try('read_images',{dirall,1:background_n*channels},mapped_molecule,handles,[]);
        for i_color=1:channels
            Temp_background=read_result(:,:,i_color:channels:size(read_result,3));
            if (i_color==2) && (channels==3)
                for i6=1:round(image_size_x/8):image_size_x
                    for i7=1:round(image_size_y/8):image_size_y
                        Temp_background1=Temp_background(i6:i6+round(image_size_x/8)-1,i7:i7+round(image_size_y/8)-1,:);
                        Temp_background1(Temp_background1>(3.*std(double(Temp_background1(:)))+ mean(Temp_background1(:))))=0;
                        Temp_background(i6:i6+round(image_size_x/8)-1,i7:i7+round(image_size_y/8)-1,:)=Temp_background1;
                    end
                end
            else
                Temp_background(Temp_background>(3.*std(double(Temp_background(:)))+ mean(Temp_background(:))))=0;
            end
            backgrounds(:,:,i_color)=sum(uint32(Temp_background),3);
            Temp_background(Temp_background~=0)=1;
            count_background(:,:,i_color)=sum(Temp_background,3);
        end
        if channels==2
            count_background(:,:,3)=uint16(ones(image_size_x,image_size_y))*settings.background_sub_max_images;
        end
        additional_required_reads=min(settings.background_sub_max_images,background_n*(background_n/min(count_background(:))-1));

        read_result=read_write_try('read_images',{dirall,background_n*channels+1:background_n*channels+additional_required_reads*channels},mapped_molecule,handles,[]);
        for i_color=1:channels
            Temp_background=read_result(:,:,i_color:3:size(read_result,3));
            if (i_color==2) && (channels==3)
                for i6=1:round(image_size_x/8):image_size_x
                    for i7=1:round(image_size_y/8):image_size_y
                        Temp_background1=Temp_background(i6:i6+round(image_size_x/8)-1,i7:i7+round(image_size_y/8)-1,:);
                        Temp_background1(Temp_background1>(3.*std(double(Temp_background1(:)))+ mean(Temp_background1(:))))=0;
                        Temp_background(i6:i6+round(image_size_x/8)-1,i7:i7+round(image_size_y/8)-1,:)=Temp_background1;
                    end
                end
            else
                Temp_background(Temp_background>(3.*std(double(Temp_background(:)))+ mean(Temp_background(:))))=0;
            end
            backgrounds(:,:,i_color)=backgrounds(:,:,i_color)+uint32(sum(Temp_background,3));
            Temp_background(Temp_background~=0)=1;
            count_background(:,:,i_color)=count_background(:,:,i_color)+uint16(sum(Temp_background,3));
            if sum(sum(count_background(:,:,i_color)==0))==0
                backgrounds(:,:,i_color)=backgrounds(:,:,i_color)./uint32(count_background(:,:,i_color));
            else
                backgrounds(:,:,i_color)=zeros(image_size_x,image_size_y);
            end        
        end
        backgrounds=uint16(backgrounds);
        Bbackground=backgrounds(:,:,1);
        if channels==3
            Gbackground=backgrounds(:,:,2);
            Rbackground=backgrounds(:,:,3);
        else if channels==2
                Gbackground=backgrounds(:,:,3);
                Rbackground=backgrounds(:,:,2);
            end
        end
        
        input1(:,:,1)=Bbackground;
        input1(:,:,2)=Gbackground;
        input1(:,:,3)=Rbackground;        
        read_write_try('write_backgrounds',input1,mapped_molecule,handles,settings) 
        clear input1;
    end
 
    rows=double(mapped_molecule{4}):double(mapped_molecule{5});
    rows_to_images=c(rows);
    [sorted_images,sorted_ind]=sort(reshape(    (rows_to_images'*ones(1,channels)+ones(length(rows_to_images),1)*(0:(channels-1)))'  ,[],1));
    read_result=read_write_try('read_images',{dirall,sorted_images},mapped_molecule,handles,[]);
    read_result=read_result(:,:,sorted_ind);    
    counter=0;
    L=[];
    Lbackground=[];
    for row=double(mapped_molecule{4}):double(mapped_molecule{5})
        B=read_result(:,:,counter+1);
        if channels==3
            G=read_result(:,:,counter+2);
            R=read_result(:,:,counter+3);
        else if channels==2
                G=uint16(zeros(image_size_x,image_size_y));
                R=read_result(:,:,counter+2);
            end
        end
        counter=counter+channels;
        L1=uint16([]);
        L1(:,:,1)=R; L1(:,:,2)=G; L1(:,:,3)=B;
        L1background(:,:,1)=Rbackground; L1background(:,:,2)=Gbackground; L1background(:,:,3)=Bbackground;
        L=color_shift_flip_stitch_and_cat(L1,L,scales,shifts,row);
        Lbackground=color_shift_flip_stitch_and_cat(L1background,Lbackground,scales,shifts,row);
    end        
    Loriginal=imrotate(L,tilt_angle.*180./pi,'bilinear','crop');
    Lbackground=imrotate(Lbackground,tilt_angle.*180./pi,'bilinear','crop');
    L=background_subtraction(Loriginal,Lbackground);
    
    ystart=max([1,1+round(double(mapped_molecule{7})-extra_length*double(mapped_molecule{10}))]);
    yend=min([size(L,1),1+round((double(mapped_molecule{5})-double(mapped_molecule{4}))*image_size_y+double(mapped_molecule{9})+sum(shifts(:,2))+extra_length*double(mapped_molecule{10}))]);
    detected_x_end=(double(mapped_molecule{5})-double(mapped_molecule{4}))*image_size_y+double(mapped_molecule{9})+sum(shifts(:,2));
    detected_x=1+(double(mapped_molecule{7}):detected_x_end);
    xmid=1+round(mean([double(mapped_molecule{6}),double(mapped_molecule{8})]));
    right_lim=min([size(L,2),xmid+xedges]);    
    left_lim=max([1,xmid-xedges]);
    molecule_image=L(:,left_lim:right_lim,:);
    molecule_image_background=Lbackground(:,left_lim:right_lim,:);
    molecule_image_original=Loriginal(:,left_lim:right_lim,:);    
    sum_range=max([1,xmid-sum_edges]):min([size(L,2),xmid+sum_edges]);
    Rstrand=sum(L(:,sum_range,1),2);
    Gstrand=sum(L(:,sum_range,2),2);
    molecule_image=molecule_image(ystart:yend,:,:);
    molecule_image_background=molecule_image_background(ystart:yend,:,:);
    molecule_image_original=molecule_image_original(ystart:yend,:,:);
    Rstrand=interp1(1:length(Rstrand),Rstrand,detected_x,'linear','extrap')';
    Gstrand=interp1(1:length(Gstrand),Gstrand,detected_x,'linear','extrap')';
    
    %calc_mol_length=((double(mapped_molecule{5})-double(mapped_molecule{4}))*image_size_y+double(mapped_molecule{9})-double(mapped_molecule{7}))+sum(shifts(:,2));
    if double(mapped_molecule{7})==0 && double(mapped_molecule{4})==1
        Rstrand=[zeros(24,1);Rstrand];
        Gstrand=[zeros(24,1);Gstrand];
        %calc_mol_length=calc_mol_length+24;
    end
    if double(mapped_molecule{9})==image_size_y-1 && double(mapped_molecule{5})==12
        Rstrand=[Rstrand;zeros(23,1)];
        Gstrand=[Gstrand;zeros(23,1)];
        %calc_mol_length=calc_mol_length+23;
    end
    %if abs(calc_mol_length-double(mapped_molecule{10}))>settings.length_calc_condition || isempty(Rstrand)
    %    Rstrand=[];
    %end
    
    if settings.write_images             
        input1{1}(:,:,:,1)=molecule_image;
        input1{1}(:,:,:,2)=molecule_image_original;
        input1{1}(:,:,:,3)=molecule_image_background;
        input1{2}=L;
        read_write_try('write_initial_images',input1,mapped_molecule,handles,settings);
    end
    %read_write_try('write_mol_length_diff',abs(calc_mol_length-double(mapped_molecule{10})),mapped_molecule,handles,settings);

else
    if isempty(dirall_pre_image_save1)
        read_result=read_write_try('load_pre_images',{dirall_pre_image_save3,dirall_pre_image_save2},mapped_molecule,handles,settings);
        molecule_image_background=read_result(:,:,:,1);
        molecule_image_original=read_result(:,:,:,2);
        molecule_image=background_subtraction(molecule_image_original,molecule_image_background);
        if settings.write_images
            read_write_try('save_image',molecule_image,mapped_molecule,handles,settings);
        end
    else
        molecule_image=read_write_try('load_final_image',dirall_pre_image_save1,mapped_molecule,handles,settings);
    end
    stitched_images_length=round((double(mapped_molecule{5})-double(mapped_molecule{4})+1)*image_size_y+sum(shifts(:,2)));
    ystart=max([1,1+round(double(mapped_molecule{7})-extra_length*double(mapped_molecule{10}))]);
    yend=min([stitched_images_length,1+round((double(mapped_molecule{5})-double(mapped_molecule{4}))*image_size_y+double(mapped_molecule{9})+sum(shifts(:,2))+extra_length*double(mapped_molecule{10}))]);
    detected_x_end=1+round((double(mapped_molecule{5})-double(mapped_molecule{4}))*image_size_y+double(mapped_molecule{9})+sum(shifts(:,2)));
    xmid=round(mean([double(mapped_molecule{6}),double(mapped_molecule{8})]));
    if xmid<xedges+1
        sum_range=size(molecule_image,2)-xedges-sum_edges:size(molecule_image,2)-xedges+sum_edges;
    else
        sum_range=xedges+1-sum_edges:min([xedges+1+sum_edges,size(molecule_image,2)]);
    end
    detected_x=(1+double(mapped_molecule{7}))-ystart:length(molecule_image)-(yend-detected_x_end);
    Rstrand=sum(molecule_image(:,sum_range,1),2);
    Gstrand=sum(molecule_image(:,sum_range,2),2);
    Rstrand=interp1(1:length(Rstrand),Rstrand,detected_x,'linear','extrap')';
    Gstrand=interp1(1:length(Gstrand),Gstrand,detected_x,'linear','extrap')';
end