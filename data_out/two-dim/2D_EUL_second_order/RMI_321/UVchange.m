clear;
UV = importdata('UV.tec',' ');
N = 51842;
Nx = 160;
Ny = 321;
mul = 8;
E = Nx * Ny;
U = UV.data((N*2+1):(N*2+E));
V = UV.data((N*2+E+1):(N*2+E*2));
U_m = reshape(U,Nx,Ny)';
V_m = reshape(V,Nx,Ny)';
Nx_p = Nx;
Ny_p = Ny;
for i=fliplr(1:Ny)
        if(mod(i,mul)~=0)
        U_m(i,:)=[];
        V_m(i,:)=[];
        Ny_p=Ny_p-1;
        end
end
for j=fliplr(1:Nx)
        if(mod(j,mul)~=0)
        U_m(:,j)=[];
        V_m(:,j)=[];    
        Nx_p=Nx_p-1;
        end
end
[X,Y] = meshgrid(linspace(0,1,Nx_p),linspace(0,2,Ny_p));
H=quiver(X,Y,U_m,V_m);
xlabel('x');
ylabel('y');
title('U-V');
axis([0,1,0,2]);
saveas(H,'UV','png')