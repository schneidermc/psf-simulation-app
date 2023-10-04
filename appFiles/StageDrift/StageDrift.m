classdef StageDrift

    % Subclasses:
    % - NoStageDrift
    % - BrownianMotion
    % - Oscillation
    % - LinearDrift
    
    properties
        motion Length
        
    end

    properties (Dependent)
        nSteps
    end

    methods
        function nSteps = get.nSteps(obj)
            nSteps = size(double(obj.motion),1);
        end

        function obj = StageDrift(motion)
            if nargin > 0
                obj.motion = motion;
            end
        end

        function obj = reduceSamplingRate(obj, nStepsSubsample)
            if mod(obj.nSteps,nStepsSubsample)~=0
                error('Number of subsampling steps must be a factor of nSteps')
            end
            subsamplingRate = obj.nSteps/nStepsSubsample;
            % lateral and axial motion are specified 
            if size(obj.motion,2)==3
                m = reshape( mean(reshape(obj.motion.inNanometer, subsamplingRate, 3, [])), obj.nSteps/subsamplingRate, [] );
                obj.motion = Length(m,'nm');
            % only lateral motion is specified 
            elseif size(obj.motion,2)==2 
                m = reshape( mean(reshape(obj.motion.inNanometer, subsamplingRate, 2, [])), obj.nSteps/subsamplingRate, [] );
                obj.motion = Length(m,'nm');
            end

  
        end
        
        % cut off motion after index cutOffPoint 
        function obj = cutOff(obj,cutOffPoint)
            cutOffPoint = floor(cutOffPoint); 
            m = obj.motion.inNanometer; 
            obj.motion = Length(m(1:cutOffPoint,:),'nm'); 
        end

        % add noise to the motion
        function obj = addNoise(obj,noiseLevel) 
            m = obj.motion.inNanometer; 
            obj.motion = Length(m + noiseLevel*randn(size(m)), 'nm'); 

        end 

        % delete axial drift component 
        function obj = deleteAxialDrift(obj)
            m = obj.motion.inNanometer; 
            obj.motion = Length(m(:,1:2),'nm'); 
        end

        % function obj = get2Ddrift(obj)



        %% Plot
        function plot(obj)
            traj = obj.motion.inNanometer;
            lineColor = lines(1);
            hold on
            plot(traj(:,1),traj(:,2),'.-','Color',lineColor)
            hStart = plot(traj(1,1),traj(1,2),'.','MarkerEdgeColor','red','MarkerFaceColor','red','MarkerSize',35);
            hEnd = plot(traj(end,1), traj(end,2),'.','Color','black','MarkerSize',35);
            axis equal
            xlabel(strcat('x [',obj.motion.getUnit,']'))
            ylabel(strcat('y [',obj.motion.getUnit,']'))
            title('Stage drift')
            legend([hStart,hEnd],{'Start','End'},'Location','northeastoutside')
            legend boxoff
        end
    end
end