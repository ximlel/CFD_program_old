line=20;
column=800;
mid=10;

L_x=1;
L_y=1;


gamma=1.4;

rho_R=1.4
rho_L=1
p_m=1
u_m=6


rho=zeros(column,1);
fid = fopen('RHO.txt','wt');
for j=1:(mid-1)
    for i=1:column
        rho(i)=rho_L;
    end
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');
end
for i=1:column
    rho(i)=rho_L+0.5*0.5*0.02*(rho_R-rho_L);
end
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
for i=1:column
    rho(i)=rho_R+0.5*0.5*0.002*(rho_L-rho_R);
end
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
for j=(mid+2):line
    for i=1:column
        rho(i)=rho_R;
    end
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');
end
fclose(fid);


u=u_m*ones(column,1);
fid = fopen('U.txt','wt');
for j=1:line
    fprintf(fid,'%12.10f\t',u);
    fprintf(fid,'\n');
end
fclose(fid);


v=zeros(column,1);
fid = fopen('V.txt','wt');
for j=1:line
    fprintf(fid,'%12.10f\t',v);
    fprintf(fid,'\n');
end
fclose(fid);


p=p_m*ones(column,1);
fid = fopen('P.txt','wt');
for j=1:line
    fprintf(fid,'%12.10f\t',p);
    fprintf(fid,'\n');
end
fclose(fid);


eps=1e-9;
t_all=100;
step=5000;

fid = fopen('config.txt','wt');
fprintf(fid,'%g\t',gamma);
fprintf(fid,'%g\t',t_all);
fprintf(fid,'%g\t',L_x);
fprintf(fid,'%g\t',L_y);
fprintf(fid,'%g\t',eps);
fprintf(fid,'%i\t',step);
fclose(fid);
