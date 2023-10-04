function [E_out,frac]=czt2(E_in,uk,ux_des,Nx)
%2D chirp z-trafo (Tong, Raoqiong, and Robert W. Cox. "Rotation of NMR images using the 2D chirp-z transform." Magnetic resonance in medicine 41.2 (1999): 253-256.
%see also OneNote (AJ)
%the normalization of E_out is consistent with the MATLAB FFT normalization
%E_in...input field
%uk...length unit in E_in
%ux_des...desired object space resolution
%Nx...size of field in pixel
%frac...fraction of energy contained in E_out in relation to the input
%field E_in
%clear all; 

% %----test parameters--------
%     Nk=64;
%     Nx=64;
%     lambda_0=532e-9;
%     RI=1;
%     NA=0.9;
%     kt_max=2*pi/lambda_0*NA;
%     uk=kt_max/(Nk/2); %k-space unit
% 
%     ux_des=100e-9;
%     [Kx,Ky,Kr,pugpil,defocus,~,uk,ux,v]=fun_usualsuspects(Nk,lambda_0,ux_des,NA,RI);
%     E_in=double(pupil);
% %---------------------------

Nk=size(E_in,1); %size of input field (k-space)
N_pad=2^nextpow2(Nx+Nk-1); %padded size required for convolution is Nx+Nk-1; here we pad more in order to have a gridsize that is a power of 2 (=faster)
N_pad=Nx+Nk-0;

x=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5))';
y=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5));

%E_in=R.^2<100^2;

%-----CZT 2D--------------------------------
    ux_fft=2*pi/(Nk*uk); %fft resolution without padding
    r=ux_des/ux_fft; %required r-factor to meet the desired resolution at the given grid size N
    alpha=r*1/Nk;
    
%   [X,Y]=ndgrid(x,y);
%   kernel2=exp(-1i*alpha*pi*(X.^2+Y.^2)); %convolution kernel
    kernel=exp(-1i*alpha*pi*x.^2)*exp(-1i*alpha*pi*y.^2); %create quadratic convolution phase kernel; faster method
    %kernel=kernel/(4*Nx.^2);
    f1=embed(E_in,[N_pad,N_pad,size(E_in,3)],0).*conj(kernel);

    %convolving f1 and f2
    %disp('czt2 time=');
    %tic
    hugo=conj(kernel).*fftshift(ifft2(fft2(ifftshift(f1)).*fft2(ifftshift(kernel))))*alpha; %multiplication with alpha leads to energy-conservation
%    E_out=embed(N_pad*conj(kernel).*hugo,Nx,0); 
    E_out=embed(hugo,[Nx,Nx,size(hugo,3)],0); 
    frac=sum(sum(abs(E_out).^2))/sum(sum(abs(E_in).^2)); %fraction of incident energy contained in E_out
    
    %imagesc(abs(hugo)); pause; 
    %toc
    
    %provisional normalization: intensities are preserved
    %E_out=E_out/sum(E_out(:).*conj(E_out(:))).*(sum(E_in(:).*conj(E_in(:))));
    
%     figure(1);
%     imagesc(abs(E_out)); axis equal; axis tight; title('CZT2'); colorbar;


% %----------comparing with FFT2--------------
%     N_pad2=round(2*pi/(ux_des*uk)); %padding required to achieve ux_des
%     disp('FFT2 time=');
%     tmp=ifftshift(embed(E_in,N_pad2,0));
%     tic
%     E=fft2(tmp);%/N_pad2; %intensity normalization
%     toc
%     figure(2);
%     E_out2=embed(fftshift(abs(E).^1),Nx,0);
%     imagesc(E_out2); axis equal; axis tight; title('FFT2'); colorbar;
%     ux_fft2=2*pi/((size(E,1))*uk);
% % 
%energy-check:
% sum(sum(abs(hugo).^2))
% sum(sum(abs(E_in).^2))