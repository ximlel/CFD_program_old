line=50;
column=100;
center=25;

drho=0;
du=0.01;
dv=0;
dp=0;

rho_c=1;
u_c=6;
v_c=0;
p_c=1;
rho_t=rho_c-drho;
u_t=u_c-du;
v_t=v_c-dv;
p_t=p_c-dp;
rho_b=rho_c+drho;
u_b=u_c+du;
v_b=v_c+dv;
p_b=p_c+dp;


fid = fopen('RHO.txt','wt');
rho=rho_c*ones(column,1);
for j=1:(center-2)
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
end
rho=rho_t*ones(column,1);
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
rho=rho_b*ones(column,1);
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
rho=rho_t*ones(column,1);
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
rho=rho_c*ones(column,1);
for j=(center+2):line
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.txt','wt');
u=u_c*ones(column,1);
for j=1:(center-2)
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
end
u=u_t*ones(column,1);
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
u=u_b*ones(column,1);
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
u=u_t*ones(column,1);
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
u=u_c*ones(column,1);
for j=(center+2):line
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.txt','wt');
v=v_c*ones(column,1);
for j=1:(center-2)
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
end
v=v_t*ones(column,1);
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
v=v_b*ones(column,1);
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
v=v_t*ones(column,1);
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
v=v_c*ones(column,1);
for j=(center+2):line
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.txt','wt');
p=p_c*ones(column,1);
for j=1:(center-2)
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
end
p=p_t*ones(column,1);
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
p=p_b*ones(column,1);
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
p=p_t*ones(column,1);
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
p=p_c*ones(column,1);
for j=(center+2):line
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
