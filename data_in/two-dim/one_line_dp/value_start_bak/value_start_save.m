line=50;
column=800;
center=25;

drho=0.001;
du=0.001;
dv=0.001;
dp=0.00001;

rho_c=1.4;
u_c=6;
v_c=0;
p_c=1;


fid = fopen('RHO.txt','wt');
rho=rho_c*ones(column,1)+drho*(rand(column,1)-0.5);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fprintf(fid,'%12.10f\t',rho_c*ones(column,1));
fprintf(fid,'\n');
for j=(center+1):line
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.txt','wt');
u=u_c*ones(column,1)+du*(rand(column,1)-0.5);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fprintf(fid,'%12.10f\t',u_c*ones(column,1));
fprintf(fid,'\n');
for j=(center+1):line
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.txt','wt');
v=v_c*ones(column,1)+dv*(rand(column,1)-0.5);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fprintf(fid,'%12.10f\t',v_c*ones(column,1));
fprintf(fid,'\n');
v=2*v_c*ones(column,1)-v;
for j=(center+1):line
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.txt','wt');
p=p_c*ones(column,1)+dp*(rand(column,1)-0.5);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fprintf(fid,'%12.10f\t',p_c*ones(column,1));
fprintf(fid,'\n');
for j=(center+1):line
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
