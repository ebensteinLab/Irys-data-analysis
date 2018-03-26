function L=red_blue_shift_uint16(L,scales)

%scales=scales./scales(1);
scales(3)=scales(3)/scales(1);
scales(2)=scales(2)/scales(1);
extra_cols=round(512.*(1-scales))./2;

for i=1:length(extra_cols)
    R=L(:,:,i);
    curr_extra_cols=extra_cols(i);
    if curr_extra_cols>=0
        if curr_extra_cols~=round(curr_extra_cols)
            R=([zeros(1,size(R,1)+1);[zeros(size(R,1),1),R]]+[[R,zeros(size(R,1),1)];zeros(1,size(R,1)+1)])./2;
            curr_extra_cols=floor(curr_extra_cols);
            %curr_extra_cols=floor(curr_extra_cols*2)/2;
        end
        Rtemp=[zeros(size(R,1),curr_extra_cols),R,zeros(size(R,1),curr_extra_cols)];
        %Rtemp=[R,zeros(size(R,1),2*curr_extra_cols)];
        %Rtemp=[zeros(size(R,1),2*curr_extra_cols),R];
        Rtemp=[zeros(curr_extra_cols,size(Rtemp,2));Rtemp;zeros(curr_extra_cols,size(Rtemp,2))];        
    else
        if curr_extra_cols~=round(curr_extra_cols)
            R=(R(2:end,2:end)+R(1:end-1,1:end-1))./2;
            curr_extra_cols=ceil(curr_extra_cols);
        end
        Rtemp=R(-curr_extra_cols+1:end+curr_extra_cols,-curr_extra_cols+1:end+curr_extra_cols);
    end
    Rtemp=imresize(Rtemp,[size(L,1),size(L,2)]);
    Rtemp(Rtemp>16383)=16383;
    Rtemp(Rtemp<0)=0;
    L(:,:,i)=Rtemp;
end