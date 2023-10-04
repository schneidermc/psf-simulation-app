classdef DoubleHelix < PhaseMask
    
    methods
        function obj = DoubleHelix(gridSize, N, d, R)
            if nargin < 1
                gridSize = 129;
            end
            if nargin < 2
                N = 9;
                d = 0.66;
                R = 1;
            end
            m = (N-1)/2;
            rVortices = (d*R) * (-m:m)';
            thetaVortices =  zeros(N,1);
            
            [theta,r] = PhaseMask.createGrid(gridSize);
            
            D = repmat(r.*exp(1i*theta),1,1,N) - reshape(rVortices.*exp(1i*thetaVortices),1,1,[]);
            M = prod(D,3);
            M = angle(M);
            M = rot90(M,1);

            obj.mask = mod(M+pi/2, 2*pi);
            obj = obj.cutAperture();
        end
    end
end