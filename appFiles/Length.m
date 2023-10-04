classdef Length < double & matlab.mixin.CustomDisplay

    properties (Access=protected)
        unit
    end

    methods
        function obj = Length(value, unit)
            obj = obj@double(value * LengthUnit.(unit))
            obj.unit = unit;
        end
        
        function unit = getUnit(obj)
            unit = obj.unit;
        end

        function val = convertUnit(obj, targetUnit)
            val = double(obj) / LengthUnit.(targetUnit);
        end

        function val = inBaseSI(obj)
            val = inMeter(obj); % base SI = meter [m]
        end

        function val = inMeter(obj)
            val = double(obj);
        end

        function val = inNanometer(obj)
            val = double(obj) / 1e-9;
        end
        
        function val = inPixels(obj, pixelSize)
            val = obj.inMeter / pixelSize.inMeter;
        end

        function objSum = plus(obj1,obj2)
            assert(strcmp(obj1.unit, obj2.unit))
            targetUnit = obj1.unit;
            sum = double(obj1)+double(obj2);
            sum = sum / LengthUnit.(targetUnit);
            objSum = Length(sum, targetUnit);
        end
    end
    
    %% custom displays 
    methods (Access=protected)
        function groups = getPropertyGroups(obj)
            info = cell(size(obj,1),1);
            val = convertUnit(obj, obj.unit);
            template = sprintf('%%8g%s ', obj.unit);
            for k = 1:size(obj,1)
                info{k} = sprintf(template, val(k,:));
            end
            groups = matlab.mixin.util.PropertyGroup(info);
        end
    end
end