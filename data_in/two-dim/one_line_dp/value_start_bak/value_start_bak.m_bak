line=50;
column=800;
center=25;

drho=0.0001*(rand-0.5);
du=0.0001*(rand-0.5);
dv=0.0001*(rand-0.5);
dp=0.0000001*(rand-0.5);

rho_c=1.4;
u_c=6;
v_c=0;
p_c=1;
rho_t=rho_c+drho;
u_t=u_c+du;
v_t=v_c+dv;
p_t=p_c+dp;
rho_b=rho_c+drho;
u_b=u_c+du;
v_b=v_c-dv;
p_b=p_c+dp;


fid = fopen('RHO.txt','wt');
rho=rho_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
rho=rho_c*ones(column,1);
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
rho=rho_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.txt','wt');
u=u_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
u=u_c*ones(column,1);
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
u=u_t*ones(column,1);

for j=(center+1):line
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.txt','wt');
v=v_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
v=v_c*ones(column,1);
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
v=v_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.txt','wt');
p=p_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
p=p_c*ones(column,1);
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
p=p_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
