line=100;
column=100;
shock=50;

gamma=1.4;


rho_R=1.4
p_R=1
u_R=10
rho_L=1
p_L=1
u_L=10

A_e = 2.5e-3
k_1 = 1.0;
k_2 = 1.0;
d_x = pi/50;
d_y = pi/50;


rho=zeros(column,1);
fid = fopen('RHO.txt','wt');
for j=1:shock
    for i=1:column
        rho(i)=rho_L;
    end
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');
end
for j=(shock+1):line
    for i=1:column
        rho(i)=rho_R;
    end
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');
end
fclose(fid);


u=zeros(column,1);
fid = fopen('U.txt','wt');
for j=1:shock
    for i=1:column
        u(i)=u_L;
    end
    fprintf(fid,'%12.10f\t',u);
    fprintf(fid,'\n');
end
for j=(shock+1):line
    for i=1:column
        u(i)=u_R;
    end
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
fid = fopen('P.txt','wt');
for j=1:shock
    for i=1:column
%        p(i)=p_L+p_L*A_e*cos(k_1*i*d_x+k_2*j*d_y);
        p(i)=p_L+p_L*A_e*(rand-0.5);
    end
    fprintf(fid,'%12.10f\t',p);
    fprintf(fid,'\n');
end
for j=(shock+1):line
    for i=1:column
%        p(i)=p_R+p_R*A_e*cos(k_1*i*d_x+k_2*j*d_y);
        p(i)=p_R+p_R*A_e*(rand-0.5);
    end
    fprintf(fid,'%12.10f\t',p);
    fprintf(fid,'\n');
end
fclose(fid);


L_x=0.01;
L_y=0.01;
eps=1e-9;
t_all=100;
step=400000;

fid = fopen('config.txt','wt');
fprintf(fid,'%g\t',gamma);
fprintf(fid,'%g\t',t_all);
fprintf(fid,'%g\t',L_x);
fprintf(fid,'%g\t',L_y);
fprintf(fid,'%g\t',eps);
fprintf(fid,'%i\t',step);
fclose(fid);
