function mapping_limits=writing_bed_and_image_files(mapping_limits,mapped_molecule,modified_x,Gstrand,Rstrand,settings)
directory=settings.directory;
i3=mapped_molecule{1};

if settings.write_bed
    Gstrand_for_bed_data=create_hmc_bed_data(mapped_molecule,Gstrand,settings.hmc_threshold,modified_x);
    if isempty(dir([directory,'bed_channel2.txt']))
        fileID3 = fopen([directory,'bed_channel2.txt'],'w');
        fprintf(fileID3,'track name="bed_channel2" description="bed_channel2" useScore=1\r\n');
        fclose(fileID3);
    end
    fileID3 = fopen([directory,'bed_channel2.txt'],'a');
    if ~isempty(Gstrand_for_bed_data)
        Gstrand_for_bed_data_str = sprintf(['chr','%d\t%16.f\t%16.f\t%d\t%d\t+\r\n'],Gstrand_for_bed_data');
        Gstrand_for_bed_data_str = regexprep(Gstrand_for_bed_data_str,'chr23','chrX');
        Gstrand_for_bed_data_str = regexprep(Gstrand_for_bed_data_str,'chr24','chrY');        
        fprintf(fileID3,Gstrand_for_bed_data_str);
    end
    fclose(fileID3);
    
    Rstrand_for_bed_data=create_hmc_bed_data(mapped_molecule,Rstrand,settings.peak_detection_threshold,modified_x);
    if isempty(dir([directory,'bed_channel1.txt']))
        fileID3 = fopen([directory,'bed_channel1.txt'],'w');
        fprintf(fileID3,'track name="bed_channel1" description="bed_channel1" useScore=1\r\n');
        fclose(fileID3);
    end
    fileID3 = fopen([directory,'bed_channel1.txt'],'a');
    if ~isempty(Rstrand_for_bed_data)
        Rstrand_for_bed_data_str = sprintf(['chr','%d\t%16.f\t%16.f\t%d\t%d\t+\r\n'],Rstrand_for_bed_data');
        Rstrand_for_bed_data_str = regexprep(Rstrand_for_bed_data_str,'chr23','chrX');
        Rstrand_for_bed_data_str = regexprep(Rstrand_for_bed_data_str,'chr24','chrY');
        fprintf(fileID3,Rstrand_for_bed_data_str);
    end
    fclose(fileID3);
    %{
    if isempty(dir([directory,'bed.gff']))
        fileID3 = fopen([directory,'bed.gff'],'w');
        fprintf(fileID3,'##gff-version 3\r\n');
        fclose(fileID3);
    end
    fileID3 = fopen([directory,'bed.gff'],'a');
    fprintf(fileID3,['chr','%d\t.\tmRNA\t%16.f\t%16.f\t%.2f\t.\t.\tID=%d\r\n'],[double(mapped_molecule{11}),modified_x(1),modified_x(end),double(mapped_molecule{17}),double(mapped_molecule{1})]);
    if ~isempty(Gstrand_for_bed_data)
        Gstrand_for_bed_data(:,4)=[];
        Gstrand_for_bed_data(:,5)=ones(size(Gstrand_for_bed_data,1),1).*double(mapped_molecule{1});
        fprintf(fileID3,['chr','%d\t.\texon\t%16.f\t%16.f\t%.2f\t.\t.\tParent=%d\r\n'],Gstrand_for_bed_data');
    end
    if ~isempty(Rstrand_for_bed_data)
        Rstrand_for_bed_data(:,4)=[];
        Rstrand_for_bed_data(:,5)=ones(size(Rstrand_for_bed_data,1),1).*double(mapped_molecule{1});
        fprintf(fileID3,['chr','%d\t.\tintron\t%16.f\t%16.f\t%.2f\t.\t.\tParent=%d\r\n'],Rstrand_for_bed_data');
    end
    fclose(fileID3);  
  %}
end

for i4=1:length(mapping_limits)+1
    if i4<=length(mapping_limits)
        if isempty(mapping_limits{i4})
            interest_zone_chr=0;
            interest_zone_bottom=0;
            interest_zone_top=0;
        else
            interest_zone_chr=mapping_limits{i4}(:,1);
            interest_zone_bottom=mapping_limits{i4}(:,2);
            interest_zone_top=mapping_limits{i4}(:,3);
        end
    else
        interest_zone_chr=0;
        interest_zone_bottom=0;
        interest_zone_top=0;
        mapping_limits{i4}=[];
    end
    if sum(mapped_molecule{11}==interest_zone_chr & ((modified_x(1)<=interest_zone_top & modified_x(end)>=interest_zone_top) | (modified_x(1)<=interest_zone_bottom & modified_x(end)>=interest_zone_bottom)  | (modified_x(1)>=interest_zone_bottom & modified_x(end)<=interest_zone_top) | (modified_x(1)<=interest_zone_bottom & modified_x(end)>=interest_zone_top)))>0
        continue;
    else
        if settings.write_graph
            if ~exist([directory,'bedgraphs_reduced\'], 'dir')
                    mkdir([directory,'bedgraphs_reduced\']);
            end
            if isempty(dir([directory,'bedgraphs_reduced\','channel1_',num2str(i4,'%d'),'.txt']))
                fileID1 = fopen([directory,'bedgraphs_reduced\','channel2_',num2str(i4,'%d'),'.txt'],'w');
                fprintf(fileID1,['track type=bedGraph name="channel2_',num2str(i4,'%d'),'" description="channel2"\r\n']);%\r\n\r');
                fclose(fileID1);
                fileID2 = fopen([directory,'bedgraphs_reduced\','channel1_',num2str(i4,'%d'),'.txt'],'w');
                fprintf(fileID2,['track type=bedGraph name="channel1_',num2str(i4,'%d'),'" description="channel1"\r\n']);%\r\n\r');
                fclose(fileID2);               
            end
            fileID1 = fopen([directory,'bedgraphs_reduced\','channel2_',num2str(i4,'%d'),'.txt'],'a');
            fileID2 = fopen([directory,'bedgraphs_reduced\','channel1_',num2str(i4,'%d'),'.txt'],'a');
            mapping_limits{i4}=[mapping_limits{i4};[mapped_molecule{11},modified_x(1),modified_x(end)]];
            midpoint_mod_x=(modified_x(1:end-1)+modified_x(2:end))./2;
            mod_x_diffs=midpoint_mod_x(2:end)'-midpoint_mod_x(1:end-1)';
            Gstrand_stretch_normalized=Gstrand(2:end-1)./mod_x_diffs.*mean(mod_x_diffs);
            Rstrand_stretch_normalized=Rstrand(2:end-1)./mod_x_diffs.*mean(mod_x_diffs);
            Gstrand_for_bed_graph=[ones(numel(modified_x)-2,1).*double(mapped_molecule{11}),round(midpoint_mod_x(1:end-1)'),round(midpoint_mod_x(2:end)'),round(Gstrand_stretch_normalized)];
            %fprintf(fileID1,['chr','%d\t%16.f\t%16.f\t%d\r\n'],Gstrand_for_bed_graph');%\r\n\n\r
            Gstrand_for_bed_graph_str = sprintf(['chr','%d\t%16.f\t%16.f\t%d\r\n'],Gstrand_for_bed_graph');
            Gstrand_for_bed_graph_str = regexprep(Gstrand_for_bed_graph_str,'chr23','chrX');
            Gstrand_for_bed_graph_str = regexprep(Gstrand_for_bed_graph_str,'chr24','chrY');
            fprintf(fileID1,Gstrand_for_bed_graph_str);
            Rstrand_for_bed_graph=[ones(numel(modified_x)-2,1).*double(mapped_molecule{11}),round(midpoint_mod_x(1:end-1)'),round(midpoint_mod_x(2:end)'),round(Rstrand_stretch_normalized)];
            %fprintf(fileID2,['chr','%d\t%16.f\t%16.f\t%d\r\n'],Rstrand_for_bed_graph');%\r\n\n\r
            Rstrand_for_bed_graph_str = sprintf(['chr','%d\t%16.f\t%16.f\t%d\r\n'],Rstrand_for_bed_graph');
            Rstrand_for_bed_graph_str = regexprep(Rstrand_for_bed_graph_str,'chr23','chrX');
            Rstrand_for_bed_graph_str = regexprep(Rstrand_for_bed_graph_str,'chr24','chrY');
            fprintf(fileID2,Rstrand_for_bed_graph_str);
            
            fclose(fileID1);
            fclose(fileID2);
        end
        runtime=mapped_molecule{19};
        runtime=runtime{:};
        ff4=strfind(runtime,'_');
        runtime=runtime(ff4(end-2)+1:end-1);
        MQR=mapped_molecule{20};
        MQR=MQR{:};
        if settings.write_images            
            %savefilename=[directory,'molecules_images\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(i3,'%d'),'_length*_not_mapped*mol.tiff'];
            savefilename=[directory,'molecules_images\*ID',num2str(i3,'%d'),'_MQR',MQR,'*.tiff'];
            dir_images=rdir(savefilename);
            [~,idx] = sort([dir_images.datenum]);
            savefilename = dir_images(idx(end)).name;
            savefilename1=savefilename(1:strfind(savefilename,'length')+5);
            savefilename1=[savefilename1,num2str(round((modified_x(end)-modified_x(1))/1000),'%dkbp'),'_bedgraph_num',num2str(i4,'%d'),'.tiff'];
            %if exist(savefilename1, 'file')
            %    delete(savefilename1);
            %end
            if ~strcmp(savefilename,savefilename1)
                movefile(savefilename,savefilename1);
            end
            %savefilename=[directory,'molecules_images\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(i3,'%d'),'_length',num2str(round((modified_x(end)-modified_x(1))/1000),'%dkbp'),'_bedgraph_num',num2str(i4,'%d'),'.tiff'];
            %mwrite(molecule_image,savefilename);
        end
        if settings.write_graph
            if ~exist([directory,'bedgraphs\'], 'dir')
                mkdir([directory,'bedgraphs\']);
            end
            savefilename=[directory,'bedgraphs\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(i3,'%d'),'_MQR',MQR,'_length',num2str(round((modified_x(end)-modified_x(1))/1000),'%dkbp'),'_bedgraph_num',num2str(i4,'%d'),'_ch2.txt'];
            fileID1 = fopen(savefilename,'w');
            fprintf(fileID1,['track type=bedGraph name="channel2_',num2str(i4,'%d'),'" description="channel2"\r\n']);
            %fprintf(fileID1,['chr','%d\t%16.f\t%16.f\t%d\r\n'],Gstrand_for_bed_graph');%\r\n\n\r
            fprintf(fileID1,Gstrand_for_bed_graph_str);
            fclose(fileID1);
            savefilename=[directory,'bedgraphs\chr',num2str(mapped_molecule{11},'%d'),'-',num2str(round(mapped_molecule{14}./1e4)./100,'%.2f'),'M-',num2str(round(mapped_molecule{15}./1e4)./100,'%.2f'),'M','_conf',num2str(round(mapped_molecule{17}),'%d'),'_run',runtime,'_scan',num2str(round(mapped_molecule{2}),'%d'),'_col',num2str(round(mapped_molecule{3}),'%d'),'_ID',num2str(i3,'%d'),'_MQR',MQR,'_length',num2str(round((modified_x(end)-modified_x(1))/1000),'%dkbp'),'_bedgraph_num',num2str(i4,'%d'),'_ch1.txt'];
            fileID2 = fopen(savefilename,'w');
            fprintf(fileID2,['track type=bedGraph name="channel1_',num2str(i4,'%d'),'" description="channel1"\r\n']);
            %fprintf(fileID2,['chr','%d\t%16.f\t%16.f\t%d\r\n'],Rstrand_for_bed_graph');%\r\n\n\r
            fprintf(fileID2,Rstrand_for_bed_graph_str);
            fclose(fileID2);
        end
        
        break;
    end
end