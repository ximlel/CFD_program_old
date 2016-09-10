line=51;
column=400;
center=26;

drho=-0.001605098804800;
du=0.005100793678606;
dv=5.524494441999758e-04;
dp=7.651804132091655e-04;

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
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
end
rho=rho_c*ones(column,1);
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
rho=rho_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.txt','wt');
u=u_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
end
u=u_c*ones(column,1);
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
u=u_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.txt','wt');
v=v_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
end
v=v_c*ones(column,1);
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
v=v_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.txt','wt');
p=p_b*ones(column,1);
for j=1:(center-1)
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
end
p=p_c*ones(column,1);
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
p=p_t*ones(column,1);
for j=(center+1):line
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
