classdef WindowZernikeAberrations_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ZernikeAberrationsUIFigure  matlab.ui.Figure
        UIAxesZernike               matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp   % Main app object
    end

    methods (Access = public)
        function initializePlot(app, zernikeAberrations)
            if nargin < 3
                zernikeAberrations = EmptyPhaseMask(129); % creates empty matrix
                zernikeAberrations = zernikeAberrations.mask;
            end

            % Plot
            imagesc(app.UIAxesZernike, zernikeAberrations);
            axis(app.UIAxesZernike, 'equal');
            axis(app.UIAxesZernike, 'tight');
            cb = colorbar(app.UIAxesZernike);
            caxis(app.UIAxesZernike, [-pi pi])
            cb.Label.String = 'Phase shift';
            cb.Ticks = (-1:0.5:1)*pi;
            cb.TickLabels = {'-\pi','','0','','\pi'};
            set(app.UIAxesZernike,'visible','off')

            app.updatePlot()
        end

        function updatePlot(app, zernikeAberrations)
            if nargin < 2
                [zernikeIndices, zernikeWeights] = app.CallingApp.readParametersZernikeAberrations();
                zernikeAberrations = ZernikeAberrations(zernikeIndices, zernikeWeights, 129);
            end
            plot(zernikeAberrations, app.UIAxesZernike);
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
        end

        % Close request function: ZernikeAberrationsUIFigure
        function ZernikeAberrationsUIFigureCloseRequest(app, event)
            app.CallingApp.ZernikeAberrationsShowplotCheckBox.Value = 0;
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ZernikeAberrationsUIFigure and hide until all components are created
            app.ZernikeAberrationsUIFigure = uifigure('Visible', 'off');
            app.ZernikeAberrationsUIFigure.Position = [1050 295 369 300];
            app.ZernikeAberrationsUIFigure.Name = 'Aberrations';
            app.ZernikeAberrationsUIFigure.CloseRequestFcn = createCallbackFcn(app, @ZernikeAberrationsUIFigureCloseRequest, true);

            % Create UIAxesZernike
            app.UIAxesZernike = uiaxes(app.ZernikeAberrationsUIFigure);
            zlabel(app.UIAxesZernike, 'Z')
            app.UIAxesZernike.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesZernike.XTick = [];
            app.UIAxesZernike.XTickLabel = '';
            app.UIAxesZernike.YTick = [];
            app.UIAxesZernike.YTickLabel = '';
            app.UIAxesZernike.FontSize = 15;
            app.UIAxesZernike.Position = [28 10 306 285];

            % Show the figure after all components are created
            app.ZernikeAberrationsUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowZernikeAberrations_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.ZernikeAberrationsUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.ZernikeAberrationsUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ZernikeAberrationsUIFigure)
        end
    end
end