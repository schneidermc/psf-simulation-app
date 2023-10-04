function [error, I_simu,P] = fun_errorcalc(Zernikes,Image_stack, alpha, apo_no, kernel, F_kernel, defocus, z_vec, E_BFP, Nk, Zernike_stack, b, norm_flag)
% Error metric calculation for vectorial phase retrieval:
% Function calculates stack of defocused bead images and compares it to a
% measured stack ("Image_stack"); output is an error metric
%
% INPUTS:
%   Zernikes    ... vector of Zernike coefficients, starting with "tip" (piston is omitted)
%   Image_stack ... measured stack of bead-images at varying values of defocus; must be normalized: sum(Image_stack(:)) = 1
%   kernel      ... convolution kernel for chirped z-trafo
%   F_kernel    ... the Fourier transform of "kernel"
%   defocus     ... phase function for 1m of defocus
%   z_vec       ... contains z-positions of images in Image_stack
%   E_BFP       ... tensor containing back focal plane fields: Ex_Px, Ex_Py, Ex_Pz, Ey_Px, Ey_Py, Ey_Pz... BFP fields
%   dia_pupil   ... size of pupil field in pixel
%   Zernike_stack   ... stack of Zernike modes to be measured
%   T_obj       ... BFP field transmission (objective apodization)
%   display     ... 'y' or 'n' if set to 'y', focal plane intensities are displayed
%   global Image_stack defocus z_vec ux uk Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz dia_pupil Zernike_stack T_obj
%   b           ... 3D-FFT of fluorescent bead; used at the very end to low-pass filter
%   norm_flag   ... if set to 1, each z-slice sum is individually normalized to 1; if set to 1, the entire stack's sum is normalized to 1
%   simulated stack
%
% OUTPUTS:
%   error   ... scalar error metric
%   I_lp    ... simulated 3D PSF (low-pass filtered by convolving with bead image)
%   a_g     ... gradient of scalar error function with respect to Zernike modes


%apo_no=2; % number of coefficients that are reserved for pupil apodization
N_img=size(Image_stack,1); 
no_images=size(Image_stack,3); % number of images in the z-stack
no_fields=size(E_BFP,3); % number of BFP-fields
Nx=size(kernel,1)-Nk; % size of simulated focal plane
%N_pad=2^nextpow2(Nx+Nk-1); % padded size required for convolution is Nx+Nk-1; here we pad more in order to have a gridsize that is a power of 2 (=faster)
N_pad=Nx+Nk-0;

% Extracting phase coefficients (Zernike-coeffs)
coefs_phase=zeros(1,1,size(Zernike_stack,3)-apo_no); %initializing
coefs_phase(1,1,(1:length(Zernikes)-apo_no))=Zernikes(1:end-apo_no); 
% Extracting amplitude apodization coefficients 
coefs_amp=zeros(1,1,apo_no);
coefs_amp(1,1,1:apo_no)=abs(Zernikes((end-apo_no+1):end));

% Handling different cases of pupil apodizations 
if apo_no==1 %no apodization 
    P=1; 
else %polynomial apodization function 
    P= sum(Zernike_stack(:,:,(end-apo_no+1):end).*repmat(coefs_amp,[Nk,Nk,1]),3);
    P(P<=0) = 0;
    %disp(P);
end
D = exp(1i*sum(Zernike_stack(:,:,1:end-apo_no).*repmat(coefs_phase,[Nk,Nk,1]),3));
conj_kernel=conj(kernel);

% Calculating stack of BFP fields, containing defocus-phases for every image plane
FF = repmat(P.*D,1,1,no_images).*exp(1i*repmat(defocus,1,1,no_images) .* permute(repmat(z_vec,Nk,1,Nk),[1,3,2]));

for q=1:no_fields
    E_simu(:,:,:,q) = czt2_fast(repmat(E_BFP(:,:,q),1,1,no_images).*FF,N_img,N_pad,alpha,kernel,F_kernel);
end

I_simu = sum(E_simu.*conj(E_simu),4);


%% Considering effects of bead-size -> low pass filtering
I_simu = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(I_simu))).*b))));

%% Normalization
if norm_flag == 1
   I_simu = I_simu ./ sum(sum(I_simu,1),2);
end

I_simu = I_simu/sum(I_simu(:));

%% Calculation of scalar error metric
error = Image_stack(:)/sum(Image_stack(:)) - I_simu(:); % for use with lsqnonlin
disp(std(error(:)));

end