A=0.142857;
K=7;







column=160;
line=321;


gamma=1.4;
Ma=1.5;
L_x=1/(column-1);
L_y=2/(line-1);


shock=round((line-1)*3/40);
center=round((line-1)/5)*2;  %center of the interface seperated two fluid.

A=A/L_y;
A_c=ceil(A);



u=0;
At=0.7575;
rho_uH=1;
rho_sH=rho_uH*(gamma+1)*Ma^2/((gamma-1)*Ma^2+2);
p_uH=1;  %p_uH=(1-1/rho_sH)*(gamma+1)/2/gamma/(Ma^2-1);  %p_L and p_uH are same.
p_sH=p_uH*(2*gamma*Ma^2-(gamma-1))/(gamma+1);
v_uH=0;
c_uH=sqrt(gamma*p_uH/rho_uH);
v_sH=(1-rho_uH/rho_sH)*c_uH*Ma+v_uH*rho_uH/rho_sH;  %v_L and v_uH are same.

rho_L=(1-At)/(1+At)*rho_uH;



index = @(y,x) y<A*cos(mod(x,column/K)*K/column*2*pi);

CC=zeros(line,column);
for j=1:(center-A_c)
CC(j,:)=ones(1,column);
end
for i=1:column
	y_1=A*cos(mod(i-1,column/K)*K/column*2*pi);
    y_2=A*cos(mod(i,column/K)*K/column*2*pi);
for j=(center-A_c+1):(center+A_c)
    if(min(y_1,y_2)>(j-center))
       CC(j,i)=1; 
    elseif(max(y_1,y_2)<(j-center-1))
       CC(j,i)=0;
    else
       CC(j,i)=quad2d(index,j-center-1,j-center,i-1,i,'AbsTol',1e-3);
    end
end
end


fid = fopen('CC.txt','wt');
for j=1:line
fprintf(fid,'%14.12f\t',CC(j,:));
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('RHO.txt','wt');
for j=1:shock
fprintf(fid,'%14.12f\t',rho_sH*ones(column,1));
fprintf(fid,'\n');
end
for j=(shock+1):line
fprintf(fid,'%14.12f\t',CC(j,:)*rho_uH+(1-CC(j,:))*rho_L);
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('U.txt','wt');
for j=1:line
fprintf(fid,'%14.12f\t',u*ones(column,1));
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('V.txt','wt');
for j=1:shock
fprintf(fid,'%14.12f\t',v_sH*ones(column,1));
fprintf(fid,'\n');
end
for j=(shock+1):line
fprintf(fid,'%14.12f\t',v_uH*ones(column,1));
fprintf(fid,'\n');
end
fclose(fid);


fid = fopen('P.txt','wt');
for j=1:shock
fprintf(fid,'%14.12f\t',p_sH*ones(column,1));
fprintf(fid,'\n');
end
for j=(shock+1):line
fprintf(fid,'%14.12f\t',p_uH*ones(column,1));
fprintf(fid,'\n');
end
fclose(fid);


eps=1e-9;
t_all=1.2;
step=25000;

fid = fopen('config.txt','wt');
fprintf(fid,'%g\t',gamma);
fprintf(fid,'%g\t',t_all);
fprintf(fid,'%g\t',L_x);
fprintf(fid,'%g\t',L_y);
fprintf(fid,'%g\t',eps);
fprintf(fid,'%i\t',step);
fclose(fid);
