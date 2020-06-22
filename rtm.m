function [img,illum]=rtm(seis,v,nbc,dx,nt,dt,s,sx,sz,gx,gz,...
    bc_top,bc_bottom,bc_left,bc_right,bc_p_nt,bc_p_nt_1)
%  IN   seis(:,:) -- seismogram,    v(:,:) -- velocity
%       nbc       -- grid number of boundary
%       dx        -- grid intervel, nt          -- number of sample
%       dt        -- time interval, s(:)        -- wavelet
%       sx,sz     -- src position,  gx(:),gz(:) -- rec position
%       bc_top(:,:,:)      -- array for storing top boundary condition
%       bc_bottom(:,:,:)   -- array for storing bottom boundary condition
%       bc_left(:,:,:)     -- array for storing left boundary condition
%       bc_right(:,:,:)    -- array for storing right boundary condition
%       bc_p_nt(:,:)       -- array for storing the last time step wavefield
%       bc_p_nt_1(:,:)     -- array for storing the second last time step wavefield
%  OUT  img(:,:)           -- Output reverse time migration image
%       illum(:,:)         --        illumination compensation

[nz,nx]=size(v); ng=numel(gx); img=zeros(nz,nx);
c=getC;
% setup ABC and temperary variables
v=padvel(v,nbc);
abc=Coef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*c(1)*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);

[isx,isz,igx,igz]=index(sx,sz,gx,gz,dx,nbc);

bp1=bc_p_nt_1;
bp0=bc_p_nt;
q0=zeros(size(v)); q1=zeros(size(v));
if nargout==2
    illum=bp1(nbc+1:nbc+nz,nbc+1:nbc+nx).^2+bp0(nbc+1:nbc+nz,nbc+1:nbc+nx).^2;
end
% Time Loop

for it=nt-2:-1:1
    % subtrace source
    bp0(isz,isx)=bp0(isz,isx)-s(it+2)*beta_dt(isz,isx);
    % dipole source
    %p0(isz-2,isx)=p0(isz-2,isx)-wavelet(it+2)*beta_dt(isz-2,isx);
    bp=temp1.*bp1-temp2.*bp0+alpha.*...
        (c(2)*(circshift(bp1,[0,1,0])+circshift(bp1,[0,-1,0])+circshift(bp1,[1,0,0])+circshift(bp1,[-1,0,0]))...
        +c(3)*(circshift(bp1,[0,2,0])+circshift(bp1,[0,-2,0])+circshift(bp1,[2,0,0])+circshift(bp1,[-2,0,0]))...
        +c(4)*(circshift(bp1,[0,3,0])+circshift(bp1,[0,-3,0])+circshift(bp1,[3,0,0])+circshift(bp1,[-3,0,0]))...
        +c(5)*(circshift(bp1,[0,4,0])+circshift(bp1,[0,-4,0])+circshift(bp1,[4,0,0])+circshift(bp1,[-4,0,0])));
    bp=load_boundary(bp,bc_top(:,:,it),bc_bottom(:,:,it),bc_left(:,:,it),bc_right(:,:,it),nz,nx,nbc); 
    q=temp1.*q1-temp2.*q0+alpha.*...
        (c(2)*(circshift(q1,[0,1,0])+circshift(q1,[0,-1,0])+circshift(q1,[1,0,0])+circshift(q1,[-1,0,0]))...
        +c(3)*(circshift(q1,[0,2,0])+circshift(q1,[0,-2,0])+circshift(q1,[2,0,0])+circshift(q1,[-2,0,0]))...
        +c(4)*(circshift(q1,[0,3,0])+circshift(q1,[0,-3,0])+circshift(q1,[3,0,0])+circshift(q1,[-3,0,0]))...
        +c(5)*(circshift(q1,[0,4,0])+circshift(q1,[0,-4,0])+circshift(q1,[4,0,0])+circshift(q1,[-4,0,0])));
    % Add seismogram
    for ig=1:ng
        q(igz(ig),igx(ig))=q(igz(ig),igx(ig))+beta_dt(igz(ig),igx(ig))*seis(it,ig);
    end
    img=image_condition(img,bp1,q0,nz,nx,nbc);
    
    if (nargout)==2
        illum=illum+bp(nbc+1:nbc+nz,nbc+1:nbc+nx).^2;
    end
    % wf refresh
    bp0=bp1; bp1=bp;
    q0=q1; q1=q;
end

end
