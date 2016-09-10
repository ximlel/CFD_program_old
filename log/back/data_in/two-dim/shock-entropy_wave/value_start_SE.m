line=32;
column=500;
shock=450;

gamma=1.4;
A_e = 0.025;
k_1=1.0;
k_2=1.0;
L_x=pi/50;
L_y=pi/16;

rho_L=1.0;
u_L=1.5;
p_L=0.714286;
rho_R=1.862069;
u_R=0.8055556;
p_R=1.755952;

rho=zeros(column,1);

fid = fopen('RHO.txt','wt');
for j=1:line
    for i=1:shock
        rho(i)=rho_L+rho_L*A_e*cos(k_1*i*L_x+k_2*(j*L_y-pi));
    end
    for i=(shock+1):column
        rho(i)=rho_R;
    end
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:shock
    u(i)=u_L;
end
for i=(shock+1):column
    u(i)=u_R;
end
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

p=zeros(column,1);
for i=1:shock
    p(i)=p_L;
end
for i=(shock+1):column
    p(i)=p_R;
end
fid = fopen('P.txt','wt');
for j=1:line
    fprintf(fid,'%12.10f\t',p);
    fprintf(fid,'\n');
end
fclose(fid);

eps=1e-9;
t_all=15;
step=5000;

fid = fopen('config.txt','wt');
fprintf(fid,'%g\t',gamma);
fprintf(fid,'%g\t',t_all);
fprintf(fid,'%g\t',L_x);
fprintf(fid,'%g\t',L_y);
fprintf(fid,'%g\t',eps);
fprintf(fid,'%i\t',step);
fclose(fid);
