function CRB = calculateCRB(obj)
    % Calculate Cramer Rao bounds for the estimation of 3D position 
    % output is Length object 
    par = obj.readParameters;
    par.shotNoise = 0; % set shotnoise to 0 for CRB calculation

    psf = PSF(par); % recalculate psf 
   
    psf = psf.image+1e-10; % add small number for numerical stability  

    % we only calculate the CRB from one z slice 
    if numel(par.defocus)>1
        nSteps = numel(par.defocus);
        midSlice = ceil(nSteps/2);
        defocusValues = par.defocus.inNanometer; 
        par.defocus = Length(defocusValues(midSlice), 'nm');
        psf = psf(:,:,midSlice); 
    end

    dx = .1; % shift in nm 
    % calculate displaced psf images;

    % we need to add the background noise in an extra step (without
    % applying Poissonian noise
    noise = par.backgroundNoise; 
    par.backgroundNoise = 0; 
   
    par.position = par.position + Length([dx 0 0],'nm');
    psf_x = PSF(par);
    par.position = par.position + Length([-dx dx 0],'nm');
    psf_y = PSF(par); 
    par.position = par.position + Length([0 -dx dx],'nm');
    psf_z = PSF(par);
    
    % calculate finite differences 
    Dx = (-psf+(psf_x.image+noise))./dx; 
    Dy = (-psf+(psf_y.image+noise))./dx; 
    Dz = (-psf+(psf_z.image+noise))./dx; 
    
    % bring matrices into same size
    LI1=1:(size(psf,1)-1);
    LI2=1:(size(psf,2)-1);
    
    Dx=Dx(LI1,LI2); 
    Dy=Dy(LI1,LI2);
    Dz=Dz(LI1,LI2); 
    
    mu=psf(LI1,LI2);
    
    
    FisherInfo=zeros(3,3);  
    
    % calculate Fisher information matrix entries
    FisherInfo(1,1)=squeeze(sum(sum(Dx.^2./mu,1),2));
    FisherInfo(2,2)=squeeze(sum(sum(Dy.^2./mu,1),2));
    FisherInfo(3,3)=squeeze(sum(sum(Dz.^2./mu,1),2));
    
    FisherInfo(1,2)=squeeze(sum(sum(Dx.*Dy./mu,1),2));
    FisherInfo(2,1)=FisherInfo(1,2);
    FisherInfo(1,3)=squeeze(sum(sum(Dx.*Dz./mu,1),2));
    FisherInfo(3,1)=FisherInfo(1,3); 
    FisherInfo(2,3)=squeeze(sum(sum(Dy.*Dz./mu,1),2));
    FisherInfo(3,2)=FisherInfo(2,3); 
    
    % calculating CRAMER RAO bound
    %inverting Fisher matrix and taking diagonal values
    FI_inv=inv(FisherInfo);
    CRBx=sqrt(FI_inv(1,1));
    CRBy=sqrt(FI_inv(2,2));
    CRBz=sqrt(FI_inv(3,3)); 

    CRB = Length([CRBx, CRBy, CRBz],'nm');

end