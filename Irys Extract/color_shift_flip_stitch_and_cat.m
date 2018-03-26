function L=color_shift_flip_stitch_and_cat(L1,L,scales,shifts,row)

L1=red_blue_shift_uint16(L1,scales);
L1=flipdim(L1 ,1); L1=flipdim(L1 ,2);

%savefilename=['e:\data\rani\temp\',num2str(row,'%d'),'.tiff'];
%imwrite(L1,savefilename);

if shifts(row,1)>0
    if shifts(row,2)~=0 %if theres a shift in both axes %shifts(row,1)*
        %make the end of previous view the average of itself and the next view
        L(end+shifts(row,2)+1:end,1+shifts(row,1):end,:)=(L(end+shifts(row,2)+1:end,1+shifts(row,1):end,:)+L1(1:-shifts(row,2),1:end-shifts(row,1),:))./2;
    end
    L1=[zeros(length(L1)+shifts(row,2),shifts(row,1),3),L1(1-shifts(row,2):end,1:end-shifts(row,1),:)];
else
    if shifts(row,2)~=0
        L(end+shifts(row,2)+1:end,1:end+shifts(row,1),:)=(L(end+shifts(row,2)+1:end,1:end+shifts(row,1),:)+L1(1:-shifts(row,2),1-shifts(row,1):end,:))./2;
    end
    L1=[L1(1-shifts(row,2):end,1-shifts(row,1):end,:),zeros(length(L1)+shifts(row,2),-shifts(row,1),3)];
end

L=cat(1,L,L1);