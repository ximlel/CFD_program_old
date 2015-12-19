line=2;
column=100;

drho=-0.0776;
du= 0.1730;
dv=0
dp=-0.0080;

rho_c=1.4;
u_c=6;
v_c=0;
p_c=1;
rho_t=rho_c+drho;
u_t=u_c+du;
v_t=v_c+dv;
p_t=p_c+dp;



fid = fopen('RHO.txt','wt');
for j=1:line
rho=rho_c*ones(column,1);
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
rho=rho_t*ones(column,1);
fprintf(fid,'%14.12f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.txt','wt');
for j=1:line
u=u_c*ones(column,1);
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
u=u_t*ones(column,1);
fprintf(fid,'%14.12f\t',u);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.txt','wt');

for j=1:line
v=v_c*ones(column,1);
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
v=v_t*ones(column,1);
fprintf(fid,'%14.12f\t',v);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.txt','wt');
for j=1:line
p=p_c*ones(column,1);
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
p=p_t*ones(column,1);
fprintf(fid,'%14.12f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
