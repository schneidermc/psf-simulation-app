function [E_BFP_x,E_BFP_y,mask_UAF,pupil]=fun_dipole_imaging(N, lambda_0, NA, RI, dipole, varargin)
%by Alexander Jesacher
%calculates the BFP-fields (Ex,Ey) of an arbitrary oriented dipole in
%the focal point of an objective
%an intermediate layer between coverglass and specimen can be assumed as
%well
%includes also near-field effects (SAF emission)
%based on the paper of Axelrod, J. of Microscopy 2012
%N...grid size
%lambda_0...vacuum wavelength
%NA...numerical aperture of objective
%RI...vector of refractive indices RI=[RI_specimen, RI_intermed., RI_immoil]
%dipole...[theta phi]   polar and azimuthal anlges of dipole
%optinal arguments: (assumed as zero if not defined)
%1) d2...thickness of intermediate layer (medium 2) in meter
%2) z...distance of dipole from medium 2
%3) f...focal length of objective (m) 
%4) mu...magnitude of dipole 
%5) T...field transmission; two pages: for s and p-pol: T=[T_s, T_p]; if
%only one matrix is provided, T_s and T_p are assumed to be the same

%% test parameters
    
%     NA=1.4; %NA of objective
%     N=64; 
%     lambda_0=532e-9; %wavelength
%     RI=[1.47, 1.47, 1.518]; %refractive indices of specimen, intermediate layer and immersion medium
%     varargin{1}=0e-9; %parameter d2 (thickness of intermediate layer (medium 2)) 
%     varargin{2}=0e-9; %z (axial position of dipole above coverslip (or intermediate layer))% 
%     varargin{3}=1e-3; %focal length of objective lens; only important for quantiatively correct E-field 
%     %defining dipole
%     varargin{4}=1e-9; %dipole magnitude
%     dipole=[0 0];
%     T = ones(N,N);

%%
nVarargs = length(varargin);

 if nVarargs==1; 
     d2=varargin{1};
     z=0; 
     f=3e-3;
     mu=1e-16;
 elseif nVarargs==2; 
     d2=varargin{1};
     z=varargin{2};
     f=3e-3;
     mu=1e-16;
 elseif nVarargs==3; 
     d2=varargin{1};
     z=varargin{2};
     f=varargin{3};
     mu=1e-16;
elseif nVarargs==4; 
     d2=varargin{1};
     z=varargin{2};
     f=varargin{3};
     mu=varargin{4};
elseif nVarargs==5 
     d2=varargin{1};
     z=varargin{2};
     f=varargin{3};
     mu=varargin{4};
     T=varargin{5}; 
     if ismatrix(T) 
         T_s=T; T_p=T; 
     else
         T_s=T(:,:,1); 
         T_p=T(:,:,2); 
     end
 else
     d2=0;
     z=0; 
     f=3e-3;
     mu=1e-16;
     T_s=1; 
     T_p=1; 
end

%% calculations
    
    theta_mu=dipole(1);  %theta-angle of dipole (0=parallel to z-axis, pi/2...in x-y-plane)
    phi_mu=dipole(2);    %phi-angle of dipole
    k0=2*pi/lambda_0;
    k=k0*RI; %k-vectors in the different media
    uk=2*k0*NA/N; %k-space unit in BFP

    if length(RI)==1; RI=[RI, RI, RI]; end
    
    %coordinates in the objective pupil
    [Kx,Ky,Kr,pupil]=create_coord(N,uk,'exact');
    PHI3=atan2(Ky,Kx);
    mask_UAF=Kr<=(k0*NA*(RI(1)/RI(3))); %BPF-zone containing UAF-light

k1=k0*RI(1); %wavevectors in media 1 (typ. water),2 (intermediate layer) and 3 (immers. medium)
k2=k0*RI(2);
k3=k0*RI(3);

THETA3=asin(Kr/k3).*pupil; %angle in medium 3  
%THETA2=acos(sqrt(1-(RI(3)/RI(2)*Kr/k3).^2)); %angle in medium2
%THETA1=acos(sqrt(1-(RI(3)/RI(1)*Kr/k3).^2)); %angle in medium1

%% calculations according to paper of Axelrod, 2012
   
% Fresnel-coefs  %Eq. 3 of paper
    CTHETA1=(sqrt(1 - (Kr/k1).^2));
    CTHETA2=(sqrt(1 - (Kr/k2).^2));
    CTHETA3=(cos(THETA3));
    
    tp12=2*RI(1)*CTHETA1./(RI(1)*CTHETA2+RI(2)*CTHETA1);
    tp23=2*RI(2)*CTHETA2./(RI(2)*CTHETA3+RI(3)*CTHETA2);
    ts12=2*RI(1)*CTHETA1./(RI(1)*CTHETA1+RI(2)*CTHETA2);
    ts23=2*RI(2)*CTHETA2./(RI(2)*CTHETA2+RI(3)*CTHETA3);
    rp12=(RI(2)*CTHETA1-RI(1)*CTHETA2)./(RI(1)*CTHETA2+RI(2)*CTHETA1);
    rp23=(RI(3)*CTHETA2-RI(2)*CTHETA3)./(RI(2)*CTHETA3+RI(3)*CTHETA2);
    rs12=(RI(1)*CTHETA1-RI(2)*CTHETA2)./(RI(1)*CTHETA1+RI(2)*CTHETA2);
    rs23=(RI(2)*CTHETA2-RI(3)*CTHETA3)./(RI(2)*CTHETA2+RI(3)*CTHETA3);
    
%from them we calculate Fresnel coefs for three-layer system
    tp=(tp12.*tp23.*exp(1i*k2*d2*CTHETA2))./(1+rp12.*rp23.*exp(2i*k2*d2*CTHETA2));
    ts=(ts12.*ts23.*exp(1i*k2*d2*CTHETA2))./(1+rs12.*rs23.*exp(2i*k2*d2*CTHETA2));

%dipole projections onto directions s, p and z:
    mu_p=mu*sin(theta_mu).*cos(phi_mu-PHI3);
    mu_s=mu*sin(theta_mu).*sin(phi_mu-PHI3);
    mu_z=mu*cos(theta_mu);

%prefactor C: (the only constant where f plays a role)
    C=(k3^2*exp(1i*k3*f).*CTHETA3)/f/RI(1).*exp(-1i*k3*d2*CTHETA3).*exp(1i*k1.*CTHETA1.*z); %eq.11 in paper
    C=C/k3.^1; %NOTE: re-normalization here, to prevent large numbers
    
%field magnitudes in layer 3 (pre-objective zone), along the s,p and z-axes (see paper for axes definitions)
    E3p=C.*tp.*CTHETA3.*(mu_p./RI(3)+mu_z.*sin(THETA3)./CTHETA1);
    E3s=C.*ts.*(mu_s/RI(3)./CTHETA1);
    E3z=C.*tp.*sin(THETA3).*(mu_p./RI(3)+mu_z.*sin(THETA3)./CTHETA1);

%influence of objective: rotation of rays by their angle theta3, such that they are all parallel to the optical axis: 
    Apo_p=1./sqrt(CTHETA3).*pupil.*T_p; %Apodization of objective lens, T_p is measured field transmission 
    Apo_s=1./sqrt(CTHETA3).*pupil.*T_s; %Apodization of objective lens, T_s is measured field transmission
    E_BFP_p=(E3p.*CTHETA3+E3z.*sin(THETA3)).*Apo_p;
    E_BFP_s=E3s.*Apo_s;  %s-polarization remains unchanged by this rotation

%coordinate transform into x-and y-polarization -> fields in the back focal plane of objective
    E_BFP_x=((cos(PHI3)).*E_BFP_p-(sin(PHI3)).*E_BFP_s); %NOTE: the final conjugation seems to be necessary (compared with Nicolas simulations), but so far I don't know why!
    E_BFP_y=((sin(PHI3)).*E_BFP_p+(cos(PHI3)).*E_BFP_s); 
      
    % display BFP intensities and phases
%     figure(1);
%     colormap gray;
%     subplot(3,2,1);
%     imagesc(abs(E_BFP_x)); axis equal; axis tight; title('|E_x|'); %colorbar;
%     subplot(3,2,2);
%     imagesc(angle(E_BFP_x)); axis equal; axis tight; title('angle(Ex)'); %colorbar;
%     subplot(3,2,3);
%     imagesc(abs(E_BFP_y)); axis equal; axis tight; title('|E_y|'); %colorbar;
%     subplot(3,2,4);
%     imagesc(angle(E_BFP_y)); axis equal; axis tight; title('angle(Ey)'); %colorbar;
%     subplot(3,2,5);
%     imagesc(abs(E_BFP_y)+abs(E_BFP_x)); axis equal; axis tight; title('|E_x|+|E_y|'); %colorbar;
    
%     figure(6);  
%     v=6;
%     quiver(X(1:v:end,1:v:end),Y(1:v:end,1:v:end),E_BFP_x(1:v:end,1:v:end),E_BFP_y(1:v:end,1:v:end));
%     axis equal; 
%     title('Electric field in BFP');
%%   
