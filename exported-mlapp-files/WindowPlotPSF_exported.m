classdef WindowPlotPSF_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PSFimageUIFigure  matlab.ui.Figure
        Toolbar           matlab.ui.container.Toolbar
        PushTool          matlab.ui.container.toolbar.PushTool
        UIAxesPSF         matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp   % Main app object
    end

    properties (Access = public)
        psfImage
    end

    methods (Access = public)
        function initializePlot(app, psfImage)
            if nargin < 2
                psfImage = zeros(17);
            end
            % PSF image
            app.psfImage = psfImage;
            imagesc(app.UIAxesPSF, psfImage);
            axis(app.UIAxesPSF, 'equal');
            axis(app.UIAxesPSF, 'tight');
            cb = colorbar(app.UIAxesPSF);
            app.setContrast(app.CallingApp.SetcontrastButtonGroup.SelectedObject.Text)
            cb.Label.String = 'Intensity';
            set(app.UIAxesPSF,'visible','off')
            colormap(app.PSFimageUIFigure, app.CallingApp.ColormapDropDown.Value)

            app.CallingApp.simulateAndDisplayPSF()
        end

        function updatePlot(app, psfImage)
            % PSF image
            app.psfImage = psfImage;
            imagesc(app.UIAxesPSF,psfImage);
            app.setContrast(app.CallingApp.SetcontrastButtonGroup.SelectedObject.Text)
            colormap(app.PSFimageUIFigure, app.CallingApp.ColormapDropDown.Value)
        end

        function setContrast(app, option)
            switch option
                case 'Intensity limits'
                    upperLimit = max([max(app.psfImage(:)), 1]);
                    caxis(app.UIAxesPSF, [0 upperLimit])
                case 'Optimize contrast'
                    lowerLimit = max([min(app.psfImage(:)),0]);
                    upperLimit = max([max(app.psfImage(:)), lowerLimit + 1]);
                    caxis(app.UIAxesPSF, [lowerLimit upperLimit])
            end
        end

        function setColormap(app,map)
            colormap(app.PSFimageUIFigure, map)
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
        end

        % Close request function: PSFimageUIFigure
        function PSFimageUIFigureCloseRequest(app, event)
            app.CallingApp.ShowPsf2DCheckBox.Value = 0;
            if app.CallingApp.ShowPsf3DCheckBox.Value == false
                app.CallingApp.ShowPSFCheckBox.Value = false;
                app.CallingApp.ShowPsf2DCheckBox.Enable = "off";
                app.CallingApp.ShowPsf3DCheckBox.Enable = "off";
            end
            app.CallingApp.CalculatingLamp.Enable = "off";
            delete(app)
        end

        % Callback function: PushTool
        function PushToolClicked(app, event)
            % save image as .mat file 
            temp = app.CallingApp.readParametersPhaseMask();
            mask = temp(129); 
            psf = app.psfImage; 
            data = {mask,psf}; 

            % save file 
            k = 0; 
            filename = ['data', num2str(k), '.mat']; 
            while isfile(filename)
                k=k+1;
                filename = ['data', num2str(k), '.mat']; 
            end
            save(filename, 'data');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PSFimageUIFigure and hide until all components are created
            app.PSFimageUIFigure = uifigure('Visible', 'off');
            app.PSFimageUIFigure.Position = [650 225 460 345];
            app.PSFimageUIFigure.Name = 'PSF image';
            app.PSFimageUIFigure.CloseRequestFcn = createCallbackFcn(app, @PSFimageUIFigureCloseRequest, true);

            % Create Toolbar
            app.Toolbar = uitoolbar(app.PSFimageUIFigure);

            % Create PushTool
            app.PushTool = uipushtool(app.Toolbar);
            app.PushTool.Tooltip = {'Save as .mat'};
            app.PushTool.ClickedCallback = createCallbackFcn(app, @PushToolClicked, true);

            % Create UIAxesPSF
            app.UIAxesPSF = uiaxes(app.PSFimageUIFigure);
            app.UIAxesPSF.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesPSF.XTick = [];
            app.UIAxesPSF.XTickLabelRotation = 0;
            app.UIAxesPSF.YTick = [];
            app.UIAxesPSF.YTickLabelRotation = 0;
            app.UIAxesPSF.ZTickLabelRotation = 0;
            app.UIAxesPSF.FontSize = 15;
            app.UIAxesPSF.Position = [30 13 368 317];

            % Show the figure after all components are created
            app.PSFimageUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowPlotPSF_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.PSFimageUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.PSFimageUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PSFimageUIFigure)
        end
    end
end