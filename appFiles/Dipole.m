classdef Dipole

    properties
        inclination (1,1) {mustBeInFullRadialRange} = 0 % Inclination angle
        azimuth (1,1) {mustBeInFullRadialRange} = 0 % Azimuthal angle
    end
    
    methods
        function obj = Dipole(angleInclination, angleAzimuth)
            obj.inclination = angleInclination;
            obj.azimuth = angleAzimuth;
        end
        
        function dipoleVector = getDipoleVector(obj)
            dipoleVector = [obj.inclination, obj.azimuth];
        end
    end
end