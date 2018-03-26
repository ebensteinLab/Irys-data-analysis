function final_image=background_subtraction(original_image,image_background)

Loriginal=original_image;
Lbackground=image_background;

str3=cell(3,1);
dirall=dir('background_equations.txt');
if ~isempty(dirall)
    str2=fileread('background_equations.txt');
else
    str2=['1:I',char(10),'2:I',char(10),'3:I'];
end
str2([strfind(str2,char(10)),strfind(str2,char(13))])=[];
str2=strrep(str2,'mean(B)','mean(B(:))');
str2=strrep(str2,'mean(I)','mean(I(:))');
str2=strrep(str2,'std(B)','std(B(:))');
str2=strrep(str2,'std(I)','std(I(:))');
str2=strrep(str2,'^','.^');
str2=strrep(str2,'/','./');
str2=strrep(str2,'*','.*');
str2=strrep(str2,'..','.');
ff1=[strfind(str2,'1:'),strfind(str2,'2:'),strfind(str2,'3:'),length(str2)+1];
for i=1:length(ff1)-1
    str3{i}=str2(ff1(i)+2:ff1(i+1)-1);
end
for i=1:3
    I=double(Loriginal(:,:,i)); B=double(Lbackground(:,:,i));
    B(B<0.5.*mean(B(:)))=mean(B(:));
    try
        L(:,:,i)=uint16(eval(str3{i}));
    catch
        L(:,:,i)=Loriginal(:,:,i);
    end
end
final_image=L;