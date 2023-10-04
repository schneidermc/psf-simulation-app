%simulates a z-stack of a fluorescent bead

function [stack_simu, E_BFP] = simu_image_stack(dia_pupil, lambda_0, NA, RI_fluid, RI_imm, dz, Nx, Nz, m_focus, ux, F_bead, slice_norm)
    
    %test parameters
%     dia_pupil = 128;
%     lambda_0 = 670e-9;
%     NA = 1.49; 
%     RI_fluid = 1.33; 
%     RI_imm = 1.518;
%     dz = 200e-9; 
%     Nx = 100; 
%     Nz = 17; 
%     m_focus = 9; 
%     ux = 100e-9;

    %fun_dipole_imaging(N,lambda_0,NA,RI,dipole,varargin)
    %calculation of BFP-fields for each dipole
    [Ex_Pz,Ey_Pz,~,~] = fun_dipole_imaging(dia_pupil,lambda_0,NA,[RI_fluid RI_fluid RI_imm],[0 0]); %z-dipole
    [Ex_Px,Ey_Px,~,~] = fun_dipole_imaging(dia_pupil,lambda_0,NA,[RI_fluid RI_fluid RI_imm],[pi/2 0]); %x-dipole
    [Ex_Py,Ey_Py,~,~] = fun_dipole_imaging(dia_pupil,lambda_0,NA,[RI_fluid RI_fluid RI_imm],[pi/2 pi/2]); %y-dipole
   
    E_BFP(:,:,1)=Ex_Px;
    E_BFP(:,:,2)=Ex_Py;
    E_BFP(:,:,3)=Ex_Pz;
    E_BFP(:,:,4)=Ey_Px;
    E_BFP(:,:,5)=Ey_Py;
    E_BFP(:,:,6)=Ey_Pz;

    %----defining parameters for chirped z-trafo----
    Nk = dia_pupil; %size of input field (in k-space)

    k0 = 2*pi/lambda_0;
    uk = 2*k0*NA/dia_pupil;

    %N_pad=2^nextpow2(Nx+Nk-1); %padded size required for convolution is Nx+Nk-1; here we pad more in order to have a gridsize that is a power of 2 (=faster)
    %N_pad = Nx + Nk-0;
    %x=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5))';
    %y=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5));
    %ux_fft=2*pi/(Nk*uk); %fft resolution without padding
    %r=ux/ux_fft; %required r-factor to meet the desired resolution at the given grid size N
    %alpha=r*1/Nk;
    %kernel=exp(-1i*alpha*pi*x.^2)*exp(-1i*alpha*pi*y.^2); %create quadratic convolution phase kernel; faster method
    %F_kernel=fft2(ifftshift(kernel));
    
    [~, ~, Kr, ~] = create_coord(dia_pupil, uk, 'exact');
    Kz = sqrt((RI_imm * k0)^2 - Kr.^2);

    %calculate 3D stack
    stack_simu = zeros(Nx,Nx,Nz); 
    I_vec = zeros(Nx,Nx,6);
    for m = 1:Nz
        for q=1:size(E_BFP,3) %propagating all BFP-fields to camera 
            I_vec(:,:,q) = abs(czt2(E_BFP(:,:,q) .* exp(1i*Kz*(m-m_focus)*dz),uk,ux,Nx)).^2;
        end    
        stack_simu(:,:,m) = sum(I_vec,3);    
    end
    
    
    %----considering effects of bead-size -> low pass filtering -----
    stack_simu = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(stack_simu))) .* F_bead))));

    %---normalization---
    if slice_norm == 1
        stack_simu = stack_simu ./ sum(sum(stack_simu,1),2); %normalize each slice
    end
    stack_simu = stack_simu/sum(stack_simu(:)); %normalize full stack

