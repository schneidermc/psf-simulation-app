function F_bead = consider_bead_size(Nx, Nz, ux, dz, m_focus, dia_bead)
           
        %-----considering effects of the bead-size----
        %simulated PSF is low-pass filtered (see Hanser et. al., J. of Microsc. 2004)
        
        %Nz=size(I4,3);
        ukz=2*pi/dz/Nz; %k-space z-coordinate unit
        ukx=2*pi/ux/Nx;
        kz=((1:Nz)-m_focus)*ukz; %kz-axis
        kx=(-Nx/2:Nx/2-1)*ukx;
        ky=kx;
        [Kx3D,Ky3D,Kz3D]=ndgrid(kx,ky,kz);
        Kr3D=sqrt(Kx3D.^2+Ky3D.^2+Kz3D.^2); %3D k-space radial coordinate
        
        bead_z = 1/2*Kr3D*dia_bead;
        F_bead = 3*(sin(bead_z)./bead_z.^3-cos(bead_z)./bead_z.^2); %Fourier Transform of the  bead shape
        F_bead(isnan(F_bead)) = 1; %remove nan in the centre