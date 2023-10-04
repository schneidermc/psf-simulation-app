function [X Y R mask]=create_coord(N,u,mode)
%creates axis with N entries (if N=[Nx Ny], a 2D grid is created)
%u...spacing between two neighbouring entries
%mode...'exact' or 'FFT'
%if 'exact' is chosen: zero point is exactly in the middle, i.e. for even grid
%sizes between the two middle pixels
%if 'FFT' is chosen: zero point is the middle pixel for odd grid sizes and
%the pixel to the lower right of the exact centre for even grid sizes.
%(this is also the pixel which MATLAB takes as zero 

if length(N)==1;
    N=[N N];
end

if length(u)==1;
    u=[u u];
end
    
    if strcmp(mode,'exact')==1;
        x=linspace(-(N(1)-1)/2,(N(1)-1)/2,N(1))*u(1);
        y=linspace(-(N(2)-1)/2,(N(2)-1)/2,N(2))*u(2);
    elseif strcmp(mode,'FFT')==1;
        x=(floor(-N(1)/2+0.5):floor(N(1)/2-0.5))*u(1);
        y=(floor(-N(2)/2+0.5):floor(N(2)/2-0.5))*u(2);
    end
    

[X,Y]=meshgrid(y,x);
R=sqrt(X.^2+Y.^2);
mask=R.^2<=((min(N/2.*u))).^2;  