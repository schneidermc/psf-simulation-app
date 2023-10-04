classdef WindowTransmissionMask_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        TransmissionmaskUIFigure  matlab.ui.Figure
        UIAxesTransmissionMask    matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp   % Main app object
    end

    methods (Access = public)
        function initializePlot(app, transmissionMask)
            if nargin < 2
                transmissionMask = EmptyTransmissionMask(129);
            end
            % Transmission mask image
            imagesc(app.UIAxesTransmissionMask, transmissionMask.mask);
            axis(app.UIAxesTransmissionMask, 'equal');
            axis(app.UIAxesTransmissionMask, 'tight');
            cb = colorbar(app.UIAxesTransmissionMask);
            colormap(app.UIAxesTransmissionMask, gray)
            caxis(app.UIAxesTransmissionMask, [0 1])
            cb.Label.String = 'Transmission';
            cb.Label.FontSize = 12;
            cb.Ticks = (0:0.25:1);
            cb.TickLabels = {'0','','0.5','','1'};
            set(app.UIAxesTransmissionMask,'visible','off')

            app.updatePlot()
        end

        function updatePlot(app, transmissionMask)
            if nargin < 2
                parTransmissionMask = app.CallingApp.readParametersTransmissionMask();
                if isa(parTransmissionMask, 'Transmission')
                    transmissionMask = parTransmissionMask;
                else
                    transmissionMask = parTransmissionMask(129);
                end
            end
            % Transmission mask image
            plot(transmissionMask, app.UIAxesTransmissionMask);
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
        end

        % Close request function: TransmissionmaskUIFigure
        function TransmissionmaskUIFigureCloseRequest(app, event)
            app.CallingApp.TransmissionMaskShowplotCheckBox.Value = 0;
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create TransmissionmaskUIFigure and hide until all components are created
            app.TransmissionmaskUIFigure = uifigure('Visible', 'off');
            app.TransmissionmaskUIFigure.Position = [1050 95 369 300];
            app.TransmissionmaskUIFigure.Name = 'Transmission mask';
            app.TransmissionmaskUIFigure.CloseRequestFcn = createCallbackFcn(app, @TransmissionmaskUIFigureCloseRequest, true);

            % Create UIAxesTransmissionMask
            app.UIAxesTransmissionMask = uiaxes(app.TransmissionmaskUIFigure);
            zlabel(app.UIAxesTransmissionMask, 'Z')
            app.UIAxesTransmissionMask.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesTransmissionMask.XTick = [];
            app.UIAxesTransmissionMask.XTickLabel = '';
            app.UIAxesTransmissionMask.YTick = [];
            app.UIAxesTransmissionMask.YTickLabel = '';
            app.UIAxesTransmissionMask.FontSize = 15;
            app.UIAxesTransmissionMask.Position = [28 10 306 285];

            % Show the figure after all components are created
            app.TransmissionmaskUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowTransmissionMask_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.TransmissionmaskUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.TransmissionmaskUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.TransmissionmaskUIFigure)
        end
    end
end