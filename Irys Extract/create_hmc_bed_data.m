function Gstrand_for_bed_data=create_hmc_bed_data(mapped_molecule,Gstrand,threshold,modified_x)

%Gstrand_mod=double(Gstrand-mean(Gstrand(Gstrand<36.*(2.*sum_edges+1))));
%[maxtab,~]=peakdet(Gstrand_mod,2.*21.*(2.*sum_edges+1),modified_x);%here is factor of 2 from the noise level
[maxtab,~]=peakdet(Gstrand,threshold,modified_x);
if ~isempty(maxtab)
    Gstrand_for_bed_data=[ones(size(maxtab,1),1).*double(mapped_molecule{11}),round(maxtab(:,1))-500,round(maxtab(:,1)+500),zeros(size(maxtab,1),1),round(maxtab(:,2))];
else
    Gstrand_for_bed_data=[];
end