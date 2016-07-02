line=100;
column=200;
shock=50;

gamma=1.4;
MS=2;

rho_star=1;
p_star=1;

rho_R=2
u_R=-0.5
v_R=1.5
p_R=0.1
rho_L=3
u_L=0.5
v_L=0.5
p_L=0.1
rho_M=2.5
u_M=-0.5*0.4+0.5*0.6
v_M=1.5*0.4+0.5*0.6
p_M=0.1

rho=zeros(column,1);
fid = fopen('RHO.txt','wt');
for j=1:line
for i=1:(line+10-j)
    rho(i)=rho_L;
end
rho(line+10-j+1)=rho_M;
for i=(line+10-j+2):column
    rho(i)=rho_R;
end
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
fid = fopen('U.txt','wt');
for j=1:line
for i=1:(line+10-j)
    u(i)=u_L;
end
u(line+10-j+1)=u_M;
for i=(line+10-j+2):column
    u(i)=u_R;
end
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.txt','wt');
for j=1:line
for i=1:(line+10-j)
    v(i)=v_L;
end
v(line+10-j+1)=v_M;
for i=(line+10-j+2):column
    v(i)=v_R;
end
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
fid = fopen('P.txt','wt');
for j=1:line
for i=1:(line+10-j)
    p(i)=p_L;
end
p(line+10-j+1)=p_M;
for i=(line+10-j+2):column
    p(i)=p_R;
end
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
