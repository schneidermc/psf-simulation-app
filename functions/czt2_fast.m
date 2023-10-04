function E_out=czt2_fast(E_in,Nx,N_pad,alpha,kernel,F_kernel)
%fast version of 2D chirped z-trafo;
%important parameters such as kernel,etc. are passed to the function

        %-----CZT 2D--------------------------------
        conj_kernel=conj(kernel);
        f1=embed(E_in,[N_pad,N_pad,size(E_in,3)],0).*conj_kernel;

        %convolving f1 and f2
        coefs_phase=conj_kernel.*fftshift(ifft2(fft2(ifftshift(f1)).*F_kernel))*alpha;
        E_out=embed(coefs_phase,[Nx,Nx,size(coefs_phase,3)],0); 