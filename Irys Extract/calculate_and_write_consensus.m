function calculate_and_write_consensus(directory,averaged_molecules,bedres,chr_lengths)

%fileID1 = fopen([directory,'bedgraphsumG.txt'],'w');
%fprintf(fileID1,['track type=bedGraph name="sumhmc" description="sumhmc"\r\n']);
%fclose(fileID1);
fileID1 = fopen([directory,'Channel2_averaged.txt'],'w');
fprintf(fileID1,['track type=bedGraph name="channel2_average" description="channel2_average"\r\n']);
fclose(fileID1);
fileID1 = fopen([directory,'Channel1_averaged.txt'],'w');
fprintf(fileID1,['track type=bedGraph name="channel1_average" description="channel1_average"\r\n']);
fclose(fileID1);
fileID1 = fopen([directory,'Molecule_count.txt'],'w');
fprintf(fileID1,['track type=bedGraph name="Molecule_count" description="Molecule_count"\r\n']);
fclose(fileID1);

fileID1 = fopen([directory,'Channel1_averaged.txt'],'a');
fileID2 = fopen([directory,'Channel2_averaged.txt'],'a');
fileID3 = fopen([directory,'Molecule_count.txt'],'a');
%fileID4 = fopen([directory,'bedgraphsumG.txt'],'a');

for i=1:24
    curr_chr=averaged_molecules{i};
    modified_x=0:bedres:chr_lengths(i);
    modified_x(curr_chr(:,3)==0)=[];
    if isempty(modified_x)
        continue;
    end
    curr_chr(curr_chr(:,3)==0,:)=[];
    %strand_for_bed_graph=[ones(length(modified_x),1).*i,modified_x',modified_x'+bedres,curr_chr(:,2)];
    %fprintf(fileID4,['chr','%d\t%16.f\t%16.f\t%d\r\n'],strand_for_bed_graph');    
    curr_chr(:,1:2)=round(curr_chr(:,1:2)./[curr_chr(:,3),curr_chr(:,3)]);    
    strand_for_bed_graph=[ones(length(modified_x),1).*i,modified_x',modified_x'+bedres,curr_chr(:,1)];
    fprintf(fileID1,['chr','%d\t%16.f\t%16.f\t%d\r\n'],strand_for_bed_graph');
    strand_for_bed_graph=[ones(length(modified_x),1).*i,modified_x',modified_x'+bedres,curr_chr(:,2)];
    fprintf(fileID2,['chr','%d\t%16.f\t%16.f\t%d\r\n'],strand_for_bed_graph');
    strand_for_bed_graph=[ones(length(modified_x),1).*i,modified_x',modified_x'+bedres,curr_chr(:,3)];
    fprintf(fileID3,['chr','%d\t%16.f\t%16.f\t%d\r\n'],strand_for_bed_graph');
end
fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
%fclose(fileID4);