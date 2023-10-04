classdef ZernikePolynomials < handle

    properties (Access=private)
        polynomials % ordering = Noll indices
        mask Mask
    end

    methods (Static)
        function obj = getInstance(mask)
            persistent instance
            if isempty(instance) || ~isequal(size(instance.getMask), size(mask)) || all(instance.getMask ~= mask,'all')
                instance = ZernikePolynomials(mask);
            end
            obj = instance;
        end
    end

    methods (Access=private)
        function obj = ZernikePolynomials(mask)
            obj.mask = mask;
            [m,n] = size(mask.values);
            obj.polynomials = NaN(m,n,0);
        end
    end

    methods (Access=public)
        function polynomials = getPolynomials(obj)
            polynomials = obj.polynomials;
        end

        function mask = getMask(obj)
            mask = obj.mask.values;
        end

        function numberCalculatedPolynomials = getNumberCalculatedPolynomials(obj)
            numberCalculatedPolynomials = size(obj.polynomials, 3);
        end

        function aberration = getAberration(obj, indicesNoll, coefficients)
            assert(numel(indicesNoll) == numel(coefficients))

            maxIndexRequired = max(indicesNoll);
            numberAlreadyCalculated = obj.getNumberCalculatedPolynomials;
            if (maxIndexRequired > numberAlreadyCalculated)
                obj.polynomials = cat(3, obj.polynomials, ...
                    obj.calculatePolynomials((numberAlreadyCalculated+1):maxIndexRequired));
            end
            % Weighted sum of Zernike polynomials
            weightedPolynomials = multiplyPagewise(obj.polynomials(:,:,indicesNoll), coefficients);
            aberration = sum( weightedPolynomials,3 );
        end

        function polynomials = calculatePolynomials(obj, nollIndices)
            [n,m] = ZernikePolynomials.noll2ZernikeIndices(nollIndices);
            polynomials = ZernikePolynomials.calculateZernikePolynomials( n, m, obj.mask, 1 );
        end
    end


    methods (Static)
        function [n,m] = noll2ZernikeIndices(j)
            n = floor(sqrt(2.*j) - 1/2);
            s = mod(n,2);
            m = (1-2.*mod(j,2)) .* (2 .* floor((1+2.*j-n.*(n+1)+s) ./ 4) - s);
        end

        function polynomials = calculateZernikePolynomials(n, m, mask, doNormalize)
            %   The radial Zernike polynomials are the radial portion of the
            %   Zernike functions, which are an orthogonal basis on the unit
            %   circle.  The series representation of the radial Zernike
            %   polynomials is
            %
            %          (n-m)/2
            %            __
            %    m      \       s                                          n-2s
            %   Z(r) =  /__ (-1)  [(n-s)!/(s!((n-m)/2-s)!((n+m)/2-s)!)] * r
            %    n      s=0
            %
            %   The following table shows the first 12 polynomials.
            %
            %       n    m    Zernike polynomial    Normalization
            %       ---------------------------------------------
            %       0    0    1                        sqrt(2)
            %       1    1    r                        sqrt(4)
            %       2    0    2*r^2 - 1                sqrt(6)
            %       2    2    r^2                      sqrt(6)
            %       3    1    3*r^3 - 2*r              sqrt(8)
            %       3    3    r^3                      sqrt(8)
            %       4    0    6*r^4 - 6*r^2 + 1        sqrt(10)
            %       4    2    4*r^4 - 3*r^2            sqrt(10)
            %       4    4    r^4                      sqrt(10)
            %       5    1    10*r^5 - 12*r^3 + 3*r    sqrt(12)
            %       5    3    5*r^5 - 4*r^3            sqrt(12)
            %       5    5    r^5                      sqrt(12)
            %       ---------------------------------------------

            %% Check input
            if nargin<5
                doNormalize = 1;
            end

            assert( all(n>=abs(m)) && all(abs(m)>=0), 'Check input! It must hold that n>=|m|>=0.')
            %assert( all(0<=rho(:)) && all(rho(:)<=1), 'Radius must be between 0 and 1!')

            %% Vectorize
            sizeMask = size(mask.values);
            [rho,theta] = mask.getPolarCoordinates();
            rho = rho(:);
            theta = theta(:);
            mAbs = abs(m);

            %% Calculate polynomials

            [requiredPowersRho, powerIdx] = ZernikePolynomials.precalculateRequiredPowersOfRho(rho, n, m);

            polynomials = NaN(sizeMask(1),sizeMask(2),numel(n));
            for j = 1:numel(n)

                A = (n(j)+mAbs(j))/2;
                D = (n(j)-mAbs(j))/2;

                %% Radiual polynomials (sum from k = 0 to (n-m)/2)

                maxFactorial = max([n(j),A]);
                cumProds = cumprod([1,1:maxFactorial]);

                % (-1)^k (n-k)!
                nominator = (-1).^(0:D) .*  cumProds(n(j)-(0:D)+1);

                % k! * ((n+m)/2-k)! * ((n-m)/2-k)!
                F = cumProds((0:D)+1);
                denominator =  F .* cumProds(A-(0:D)+1) .* fliplr(F);

                % r^(n-2k)
                powers = n(j)-2*(0:D);
                powersRho = requiredPowersRho(:,powerIdx(powers+1));

                Rnm = sum(nominator./denominator .* powersRho,2);

                %% Zernike polynomials
                if m(j)==0
                    Z = Rnm;
                elseif m(j)>0 % 'even' Zernike polynomials
                    Z = Rnm .* cos(theta*mAbs(j));
                else % 'odd' Zernike polynomials
                    Z = Rnm .* sin(theta*mAbs(j));
                end

                %% Normalization
                if doNormalize
                    if m(j)==0
                        Z = Z *sqrt(n(j)+1);
                    else
                        Z = Z * sqrt(2*(n(j)+1));
                    end
                end

                polynomials(:,:,j) = reshape(Z,sizeMask).*mask.values;
            end
        end

        function [requiredPowersRho, powerIdx] = precalculateRequiredPowersOfRho(rho, n, m)
            isEven = mod(n,2);
            if all(isEven) || all(~isEven)
                requiredPowers = min(abs(m)):2:max(n);
            else
                requiredPowers = min(abs(m)):1:max(n);
            end
            if requiredPowers(1)==0 % faster version if power 0 is required
                requiredPowersRho = arrayfun(@(p) rho.^p, requiredPowers(2:end), 'UniformOutput', false);
                requiredPowersRho = cat(2, requiredPowersRho{:});
                requiredPowersRho = [ones(length(rho),1) requiredPowersRho];
            else
                requiredPowersRho = arrayfun(@(p) rho.^p, requiredPowers, 'UniformOutput',false);
                requiredPowersRho = cat(2, requiredPowersRho{:});
            end
            powerIdx = NaN(max(n)+1,1);
            powerIdx(requiredPowers+1) = 1:length(requiredPowers);
        end
    end
end