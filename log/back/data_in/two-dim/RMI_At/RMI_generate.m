for j=1:20
    A=linspace(0.1,0.9,20);
    i=A(j);
    foldername=strcat(' RMI_At_',num2str(j,'%g'));
    s=strcat('mkdir',foldername);
    system(s);
	foldername=strcat('RMI_At_',num2str(j,'%g'));
    filename=strcat(foldername,'/value_start.m');
    copyfile('value_start.m_bak',filename)
	fid = fopen(filename,'rt+');
    frewind(fid);
    fprintf(fid,'At=%g;\n',i);
    fclose(fid);
end
