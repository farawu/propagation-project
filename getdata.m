clear all; 

cLight = 3e8;
frq = 300e+06; 
lamda = cLight/frq;
X=1000;
Y=1000;
Z=6000;
dx=1;
dy=1;
dz=5;
Xmax=X./dx;
Ymax=Y./dy;  %
Zmax=Z./dz;% progation direction


er=h5read('sbopaperout.h5','/yzmapr');
ei=h5read('sbopaperout.h5','/yzmapi');
e=er+1i*ei;
Uyz=zeros(Zmax,Ymax);
for i=1:Ymax
    Uyz(:,Ymax-i+1)=e((i-1)*Zmax+1:i*Zmax);
end



er2=h5read('sbopaperout.h5','/xzmapr');
ei2=h5read('sbopaperout.h5','/xzmapi');
e2=er2+1i*ei2;
Uxz=zeros(Zmax,Xmax);
for i=1:Xmax
    Uxz(:,i)=e2((i-1)*Zmax+1:i*Zmax);
end

er3=h5read('sbopaperout.h5','/xymapr');
ei3=h5read('sbopaperout.h5','/xymapi');
e3=er3+1i*ei3;
Uxy=zeros(Ymax,Xmax);
for i=1:Xmax
    Uxy(:,i)=e3((i-1)*Ymax+1:i*Ymax);
end

pathloss_yz=zeros(Zmax,Ymax);
Zz=(1:Zmax)*dz;
pathloss_zx=meshgrid(Zz,1:Ymax)';

pathloss_xz=zeros(Zmax,Xmax);
Zz=(1:Zmax)*dz;
pathloss_xz=meshgrid(Zz,Xmax)';




 %% YZ面

dyt =Y:dy:0;
dzt = 0:dz:Z;
koko=abs(Uyz);
 koko2=-20*log10(max(koko,1e-10))+20*log10(4*pi)-30*log10(lamda)+10.*log10(pathloss_zx);%
 Uro=rot90(koko2);%rot90
 Uabs=Uro;
 figure;
image(dzt,dyt,Uabs,'CDataMapping','scaled')%
colormap(jet);
colorbar('EastOutside');
 set(gca,'YDir','normal');
title('YZ-cross section');
xlabel('Z');
ylabel('Y');


 %% XZ面

dzt = 1:dz:Z;
dxt = -X:dx:X;
koko=abs(Uxz);
koko2=-20*log10(max(koko,1e-10))+20*log10(4*pi)-30*log10(lamda)+10.*log10(pathloss_xz);%
Uro=rot90(koko2);
Uabs=Uro;
 figure;
image(dzt,dxt,Uabs,'CDataMapping','scaled')
colormap(jet);
colorbar('EastOutside'); 
title('XZ-cross section');
xlabel('Z');
ylabel('X');

 %% XY面

dyt = 0:dy:Y;
dxt = -X:dx:X;
koko=abs(Uxy);
koko2=-20*log10(max(koko,1e-10))+20*log10(4*pi)-30*log10(lamda);%+10.*log10(pathloss_xz)
Uro=fliplr(koko2);
Uabs=Uro;
 figure;
image(dxt,dyt,Uabs,'CDataMapping','scaled')
colormap(jet);
colorbar('EastOutside'); 
set(gca,'YDir','normal');
title('XY-cross section');
xlabel('Y');
ylabel('X');


