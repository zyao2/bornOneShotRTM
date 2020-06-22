function seis=forward(v,refl,nbc,dx,nt,dt,s,sx,sz,gx,gz)
%  IN   v(:,:) -- velocity,      refl(:,:)   -- reflectivity 
%       nbc    -- grid number of boundary
%       dx     -- grid intervel, nt          -- number of sample
%       dt     -- time interval, s(:)        -- wavelet
%       sx,sz  -- src position,  gx(:),gz(:) -- rec position
%       isFS   -- if is Free Surface condition
%  OUT  seis(:,:)          -- Output born modeling seismogram
%       bc_top(:,:,:)      -- array for storing top boundary condition
%       bc_bottom(:,:,:)   -- array for storing bottom boundary condition
%       bc_left(:,:,:)     -- array for storing left boundary condition
%       bc_right(:,:,:)    -- array for storing right boundary condition
%       bc_p_nt(:,:)       -- array for storing the last time step wavefield
%       bc_p_nt_1(:,:)     -- array for storing the second last time step wavefield

seis=zeros(nt,numel(gx));[nz,nx]=size(v);x = (0:nx-1)*dx; z = (0:nz-1)*dx;
ng=numel(gx); g=1:ng; t=(0:nt-1)*dt;
c=getC;
% setup ABC and temperary variables
v=padvel(v,nbc); refl=padvel(refl,nbc);
abc=Coef2D(v,nbc,dx);
alpha=(v*dt/dx).^2; kappa=abc*dt;
temp1=2+2*c(1)*alpha-kappa; temp2=1-kappa;
beta_dt = (v*dt).^2;
s=expand_source(s,nt);
[isx,isz,igx,igz]=index(sx,sz,gx,gz,dx,nbc);

p1=zeros(size(v)); p0=zeros(size(v));
q1=zeros(size(v)); q0=zeros(size(v));

% Time Looping
for it=1:nt
    p=temp1.*p1-temp2.*p0+alpha.*...
        (c(2)*(circshift(p1,[0,1,0])+circshift(p1,[0,-1,0])+circshift(p1,[1,0,0])+circshift(p1,[-1,0,0]))...
        +c(3)*(circshift(p1,[0,2,0])+circshift(p1,[0,-2,0])+circshift(p1,[2,0,0])+circshift(p1,[-2,0,0]))...
        +c(4)*(circshift(p1,[0,3,0])+circshift(p1,[0,-3,0])+circshift(p1,[3,0,0])+circshift(p1,[-3,0,0]))...
        +c(5)*(circshift(p1,[0,4,0])+circshift(p1,[0,-4,0])+circshift(p1,[4,0,0])+circshift(p1,[-4,0,0])));
    p(isz,isx) = p(isz,isx) + beta_dt(isz,isx) * s(it);
 
    q=temp1.*q1-temp2.*q0+alpha.*...
        (c(2)*(circshift(q1,[0,1,0])+circshift(q1,[0,-1,0])+circshift(q1,[1,0,0])+circshift(q1,[-1,0,0]))...
        +c(3)*(circshift(q1,[0,2,0])+circshift(q1,[0,-2,0])+circshift(q1,[2,0,0])+circshift(q1,[-2,0,0]))...
        +c(4)*(circshift(q1,[0,3,0])+circshift(q1,[0,-3,0])+circshift(q1,[3,0,0])+circshift(q1,[-3,0,0]))...
        +c(5)*(circshift(q1,[0,4,0])+circshift(q1,[0,-4,0])+circshift(q1,[4,0,0])+circshift(q1,[-4,0,0])));
    q=pertubation(p1,q,refl,beta_dt);
    for ig=1:ng
        seis(it,ig)=q(igz(ig),igx(ig));
    end

    p0=p1; p1=p;
    q0=q1; q1=q;
end
