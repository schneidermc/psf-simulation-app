classdef WindowPhaseMask_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PhasemaskUIFigure  matlab.ui.Figure
        UIAxesPhaseMask    matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp   % Main app object
    end

    methods (Access = public)
        function initializePlot(app, phaseMask)
            if nargin < 2
                phaseMask = EmptyPhaseMask(129);
            end
            % Phase mask image
            imagesc(app.UIAxesPhaseMask, phaseMask.mask);
            axis(app.UIAxesPhaseMask, 'equal');
            axis(app.UIAxesPhaseMask, 'tight');
            cb = colorbar(app.UIAxesPhaseMask);
            if isMATLABReleaseOlderThan('R2022a')
                caxis(app.UIAxesPhaseMask, [0 2*pi])
            else
                clim(app.UIAxesPhaseMask, [0 2*pi])
            end
            cb.Label.String = 'Phase shift';
            cb.Ticks = (0:0.5:2)*pi;
            cb.TickLabels = {'0','','\pi','','2\pi'};
            set(app.UIAxesPhaseMask,'visible','off')

            app.updatePlot()
        end

        function updatePlot(app, phaseMask)
            if nargin < 2
                parPhaseMask = app.CallingApp.readParametersPhaseMask();
                if isa(parPhaseMask,'PhaseMask')
                    phaseMask = parPhaseMask;
                else
                    phaseMask = parPhaseMask(129);
                end
            end
            % Phase mask image
            plot(phaseMask, app.UIAxesPhaseMask);
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
        end

        % Close request function: PhasemaskUIFigure
        function PhasemaskUIFigureCloseRequest(app, event)
            app.CallingApp.PhaseMaskShowplotCheckBox.Value = 0;
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PhasemaskUIFigure and hide until all components are created
            app.PhasemaskUIFigure = uifigure('Visible', 'off');
            app.PhasemaskUIFigure.Position = [1050 295 369 300];
            app.PhasemaskUIFigure.Name = 'Phase mask';
            app.PhasemaskUIFigure.CloseRequestFcn = createCallbackFcn(app, @PhasemaskUIFigureCloseRequest, true);

            % Create UIAxesPhaseMask
            app.UIAxesPhaseMask = uiaxes(app.PhasemaskUIFigure);
            zlabel(app.UIAxesPhaseMask, 'Z')
            app.UIAxesPhaseMask.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesPhaseMask.XTick = [];
            app.UIAxesPhaseMask.XTickLabel = '';
            app.UIAxesPhaseMask.YTick = [];
            app.UIAxesPhaseMask.YTickLabel = '';
            app.UIAxesPhaseMask.FontSize = 15;
            app.UIAxesPhaseMask.Position = [28 10 306 285];

            % Show the figure after all components are created
            app.PhasemaskUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowPhaseMask_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.PhasemaskUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.PhasemaskUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PhasemaskUIFigure)
        end
    end
end