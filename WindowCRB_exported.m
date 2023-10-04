classdef WindowCRB_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        CramerRaoboundUIFigure  matlab.ui.Figure
        UIAxesCRBx              matlab.ui.control.UIAxes
        UIAxesCRB               matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp   % Main app object
    end

    methods (Access = public)
        function updatePlot(app, psf)
            % Calculate CRB
            n = 30;
            CRBx = NaN(n,1); CRBy = NaN(n,1); CRBd = NaN(n,1);
            CRBtheta = NaN(n,1); CRBphi = NaN(n,1);
            thetaVec = zeros(n,1);
            for k=1:n
                thetaVec(k) = (k-1/2)*pi / (2*n);
                [CRBx(k), CRBy(k), CRBd(k), CRBtheta(k), CRBphi(k)] = boundsCramerRao(psf, thetaVec(k));
            end
            % Plot graphs
            plot(app.UIAxesCRB, thetaVec,180/pi*sqrt(CRBtheta), thetaVec,180/pi*sqrt(CRBphi))
            legend(app.UIAxesCRB, {'\theta', '\phi'})

            plot(app.UIAxesCRBx,thetaVec,sqrt(CRBx), thetaVec,sqrt(CRBy))
            legend(app.UIAxesCRBx, {'x', 'y'})
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
        end

        % Close request function: CramerRaoboundUIFigure
        function CramerRaoboundUIFigureCloseRequest(app, event)
            % Close window
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create CramerRaoboundUIFigure and hide until all components are created
            app.CramerRaoboundUIFigure = uifigure('Visible', 'off');
            app.CramerRaoboundUIFigure.Position = [100 80 610 252];
            app.CramerRaoboundUIFigure.Name = 'Cramer-Rao bound';
            app.CramerRaoboundUIFigure.CloseRequestFcn = createCallbackFcn(app, @CramerRaoboundUIFigureCloseRequest, true);

            % Create UIAxesCRB
            app.UIAxesCRB = uiaxes(app.CramerRaoboundUIFigure);
            title(app.UIAxesCRB, 'CRB')
            xlabel(app.UIAxesCRB, '\theta')
            ylabel(app.UIAxesCRB, 'CRB [°]')
            app.UIAxesCRB.PlotBoxAspectRatio = [1.30615164520744 1 1];
            app.UIAxesCRB.XLim = [0 1.5707963267949];
            app.UIAxesCRB.XTick = [0 0.392699081698724 0.785398163397448 1.17809724509617 1.5707963267949];
            app.UIAxesCRB.XTickLabelRotation = 0;
            app.UIAxesCRB.XTickLabel = {'0'; '\pi/8'; '\pi/4'; '3\pi/8'; '\pi/2'};
            app.UIAxesCRB.YTickLabelRotation = 0;
            app.UIAxesCRB.ZTickLabelRotation = 0;
            app.UIAxesCRB.BusyAction = 'cancel';
            app.UIAxesCRB.Position = [14 4 273 234];

            % Create UIAxesCRBx
            app.UIAxesCRBx = uiaxes(app.CramerRaoboundUIFigure);
            title(app.UIAxesCRBx, 'CRB')
            xlabel(app.UIAxesCRBx, '\theta')
            ylabel(app.UIAxesCRBx, 'CRB [nm]')
            app.UIAxesCRBx.PlotBoxAspectRatio = [1.30615164520744 1 1];
            app.UIAxesCRBx.XLim = [0 1.5707963267949];
            app.UIAxesCRBx.XTick = [0 0.392699081698724 0.785398163397448 1.17809724509617 1.5707963267949];
            app.UIAxesCRBx.XTickLabelRotation = 0;
            app.UIAxesCRBx.XTickLabel = {'0'; '\pi/8'; '\pi/4'; '3\pi/8'; '\pi/2'};
            app.UIAxesCRBx.YTickLabelRotation = 0;
            app.UIAxesCRBx.ZTickLabelRotation = 0;
            app.UIAxesCRBx.BusyAction = 'cancel';
            app.UIAxesCRBx.Position = [307 4 273 234];

            % Show the figure after all components are created
            app.CramerRaoboundUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowCRB_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.CramerRaoboundUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.CramerRaoboundUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.CramerRaoboundUIFigure)
        end
    end
end