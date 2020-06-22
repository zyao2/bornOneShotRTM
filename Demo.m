%Born modeling and adjoint
% the adjoint operation of RTM, adjoint test and migration green's function.
%velocity model
clear all; close all
nz=81;nx=201;
vel=zeros(nz,nx);
vel(1:30,:)=1000;
vel(31:60,:)=1200;
vel(61:end,:)=1500;

%source and receiver location;
dx=5;dz=5;
sx=100*dx;sz=0;
recx=(0:2:(nx-1))*dx; recz=zeros(size(recx));
x = (0:nx-1)*dx; z = (0:nz-1)*dx;

%FD parameters;
nbc=20; nt=2001; dt=0.0005;

%source wavelet;
freq=25; s=ricker(freq,dt);
%Smooth the true veolicyt to get the migration velocity model;

[vel_ss,refl_ss]=vel_smooth(vel,3,3,1);

%Plot the background velocity, reflectivity model and wavelet;

figure;set(gcf,'position',[0 0 600 300]);colormap(gray);
subplot(131);imagesc(x,z,vel_ss);colorbar;
xlabel('X (m)'); ylabel('Z (m)'); title('smooth velocity');
figure(1);subplot(132);imagesc(x,z,refl_ss);colorbar;
xlabel('X (m)'); ylabel('Z (m)'); title('reflectivity');
figure(1);subplot(133);plot((0:numel(s)-1)*dt,s);
xlabel('Time (s)'); ylabel('Amplitude');title('wavelet');
drawnow
%Run the modeling code, watch the movie of wave propagation including both 
%background wave field and pertubation wave field, then plot the seismic data;
%tic;
%seis=forward(vel_ss,refl_ss,nbc,dx,nt,dt,s,sx,sz,recx,recz);
%save seis seis
%toc;
load seis  
%Do the RTM with the generated born data;

disp('Forward Modeling to save BC')

tic; 
[~,bc_top,bc_bottom,bc_left,bc_right,bc_p_nt,bc_p_nt_1]=...
rtm_forwd(vel_ss,nbc,dx,nt,dt,s,sx,sz,recx,recz);toc;
disp(' RTM ')
img=rtm(seis,vel_ss,nbc,dx,nt,dt,s,sx,sz,recx,recz,...
bc_top,bc_bottom,bc_left,bc_right,bc_p_nt,bc_p_nt_1);
figure;set(gcf,'position',[0 0 600 300]);colormap(gray);
imagesc(x,z,img);caxis([-100 100]);
xlabel('X (m)'); ylabel('Z (m)'); title('RTM Image with Born Data');
