classdef WindowPlotPSFThreeDim_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PSFThreeDimUIFigure       matlab.ui.Figure
        Toolbar                   matlab.ui.container.Toolbar
        PushTool                  matlab.ui.container.toolbar.PushTool
        TransparencySpinner       matlab.ui.control.Spinner
        TransparencySpinnerLabel  matlab.ui.control.Label
        UpdatelightingButton      matlab.ui.control.Button
        IsoSurfaceSlider          matlab.ui.control.Slider
        IsosurfaceLabel           matlab.ui.control.Label
        UIAxesIsoSurface          matlab.ui.control.UIAxes
        UIAxesPSF                 matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp % Main app object
        isosurfacePatch
        cameraLight
        isosurfaceTransparency = 0.2
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
            %axis(app.UIAxesPSF, 'equal');
            axis(app.UIAxesPSF, 'tight');
            cb = colorbar(app.UIAxesPSF);
            %app.setContrast(app.CallingApp.SetcontrastButtonGroup.SelectedObject.Text)
            cb.Label.String = 'Intensity';
            set(app.UIAxesPSF,'visible','on')
            colormap(app.PSFThreeDimUIFigure, app.CallingApp.ColormapDropDown.Value)

            app.CallingApp.simulateAndDisplayPSF();
        end

        function updatePlot(app, psfImage)
            % PSF image
            app.psfImage = psfImage;

            % Projected side view (xz):
            projection = squeeze(sum(psfImage,1));
            imagesc(app.UIAxesPSF,projection);
            
            nPixels = app.CallingApp.nPixels;
            pixelSize = app.CallingApp.PixelsizeObjectSpaceSpinner.Value;
            xAxis = linspace(0, pixelSize*nPixels, nPixels) - pixelSize*nPixels/2;
            zAxis = app.CallingApp.parameters.defocus.inNanometer();
            
            % Show xz projection
            imagesc(app.UIAxesPSF, xAxis, zAxis, projection');
            zLimit = 1/2*app.CallingApp.StepsizeEditField.Value*(app.CallingApp.NumberstepsEditField.Value-1);
            app.UIAxesPSF.YTick = [zAxis(1), 0, zAxis(length(zAxis))];
            app.UIAxesPSF.YTickLabel = {num2str(-zLimit), num2str(0), num2str(zLimit)};
            colormap(app.UIAxesPSF, app.CallingApp.ColormapDropDown.Value)

            % Calculate isosurface 
            calculateIsoSurface(app,psfImage);
        end


        function calculateIsoSurface(app,psfImage)
            maxIntensity = max(psfImage(:)); 
            app.IsoSurfaceSlider.Limits = [0 1];
            app.IsoSurfaceSlider.MajorTicks = [0 1];
            app.IsoSurfaceSlider.MajorTickLabels = {'0', num2str(100*ceil(maxIntensity/100))};

            % Isosurface      
            dim = size(psfImage);
            x = 1:dim(1); 
            y = 1:dim(2); 
            z = 1:dim(3); 

            cmap = get(app.UIAxesPSF,'colormap');

            cla(app.UIAxesIsoSurface)

            [faces,verts] = isosurface(x,y,z,psfImage,max(1, 100*ceil(maxIntensity/100)*app.IsoSurfaceSlider.Value));
            app.isosurfacePatch = patch(app.UIAxesIsoSurface, 'Vertices',verts,'Faces',faces, 'FaceColor',cmap(max(1,ceil(256*app.IsoSurfaceSlider.Value)),:), 'EdgeColor', 'none');

            view(app.UIAxesIsoSurface, 3)
            app.UIAxesIsoSurface.XLim = [1 dim(1)];
            app.UIAxesIsoSurface.YLim = [1 dim(2)];
            app.UIAxesIsoSurface.ZLim = [1 dim(3)];

            % Lighting
            app.isosurfacePatch.FaceVertexAlphaData = app.isosurfaceTransparency;
            app.isosurfacePatch.FaceAlpha = 'flat';
            lighting(app.UIAxesIsoSurface, 'gouraud');
            material(app.UIAxesIsoSurface, 'dull');
            app.cameraLight = camlight(app.UIAxesIsoSurface, 'right');
            
            % Single xy plane:
            %imagesc(app.UIAxesPSF,psfImage(:,:,1));
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
            colormap(app.PSFThreeDimUIFigure, map)
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
        end

        % Close request function: PSFThreeDimUIFigure
        function PSFThreeDimUIFigureCloseRequest(app, event)
            app.CallingApp.ShowPsf3DCheckBox.Value = 0;
            if app.CallingApp.ShowPsf2DCheckBox.Value == false
                app.CallingApp.ShowPSFCheckBox.Value = false;
                app.CallingApp.ShowPsf2DCheckBox.Enable = "off";
                app.CallingApp.ShowPsf3DCheckBox.Enable = "off";
            end
            delete(app)
        end

        % Clicked callback: PushTool
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

        % Value changed function: IsoSurfaceSlider
        function IsoSurfaceSliderValueChanged(app, event)
            app.IsoSurfaceSlider.Value = event.Value;
            app.calculateIsoSurface(app.psfImage);
        end

        % Value changing function: IsoSurfaceSlider
        function IsoSurfaceSliderValueChanging(app, event)
            app.IsoSurfaceSlider.Value = event.Value;
            app.calculateIsoSurface(app.psfImage);
        end

        % Button pushed function: UpdatelightingButton
        function UpdatelightingButtonPushed(app, event)
            camlight(app.cameraLight, 'right');
        end

        % Value changed function: TransparencySpinner
        function TransparencySpinnerValueChanged(app, event)
            app.TransparencySpinner.Value = event.Value;
            app.isosurfaceTransparency = app.TransparencySpinner.Value;
            alpha(app.UIAxesIsoSurface, app.isosurfaceTransparency);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PSFThreeDimUIFigure and hide until all components are created
            app.PSFThreeDimUIFigure = uifigure('Visible', 'off');
            app.PSFThreeDimUIFigure.Position = [650 225 772 360];
            app.PSFThreeDimUIFigure.Name = 'PSF 3D';
            app.PSFThreeDimUIFigure.CloseRequestFcn = createCallbackFcn(app, @PSFThreeDimUIFigureCloseRequest, true);

            % Create Toolbar
            app.Toolbar = uitoolbar(app.PSFThreeDimUIFigure);

            % Create PushTool
            app.PushTool = uipushtool(app.Toolbar);
            app.PushTool.Tooltip = {'Save as .mat'};
            app.PushTool.ClickedCallback = createCallbackFcn(app, @PushToolClicked, true);

            % Create UIAxesPSF
            app.UIAxesPSF = uiaxes(app.PSFThreeDimUIFigure);
            title(app.UIAxesPSF, 'xz-Projection')
            ylabel(app.UIAxesPSF, 'z [nm]')
            app.UIAxesPSF.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesPSF.XLim = [0 1];
            app.UIAxesPSF.XTick = [];
            app.UIAxesPSF.XTickLabel = '';
            app.UIAxesPSF.YTick = [];
            app.UIAxesPSF.FontSize = 12;
            app.UIAxesPSF.Position = [48 61 269 272];

            % Create UIAxesIsoSurface
            app.UIAxesIsoSurface = uiaxes(app.PSFThreeDimUIFigure);
            title(app.UIAxesIsoSurface, 'Isosurface')
            app.UIAxesIsoSurface.XTick = [];
            app.UIAxesIsoSurface.YTick = [];
            app.UIAxesIsoSurface.ZTick = [];
            app.UIAxesIsoSurface.BusyAction = 'cancel';
            app.UIAxesIsoSurface.Position = [402 61 259 272];

            % Create IsosurfaceLabel
            app.IsosurfaceLabel = uilabel(app.PSFThreeDimUIFigure);
            app.IsosurfaceLabel.HorizontalAlignment = 'center';
            app.IsosurfaceLabel.Position = [660 277 93 23];
            app.IsosurfaceLabel.Text = 'Surface value';

            % Create IsoSurfaceSlider
            app.IsoSurfaceSlider = uislider(app.PSFThreeDimUIFigure);
            app.IsoSurfaceSlider.Limits = [0 1];
            app.IsoSurfaceSlider.MajorTicks = [0 1];
            app.IsoSurfaceSlider.MajorTickLabels = {''};
            app.IsoSurfaceSlider.Orientation = 'vertical';
            app.IsoSurfaceSlider.ValueChangedFcn = createCallbackFcn(app, @IsoSurfaceSliderValueChanged, true);
            app.IsoSurfaceSlider.ValueChangingFcn = createCallbackFcn(app, @IsoSurfaceSliderValueChanging, true);
            app.IsoSurfaceSlider.BusyAction = 'cancel';
            app.IsoSurfaceSlider.Position = [715 111 3 150];
            app.IsoSurfaceSlider.Value = 0.5;

            % Create UpdatelightingButton
            app.UpdatelightingButton = uibutton(app.PSFThreeDimUIFigure, 'push');
            app.UpdatelightingButton.ButtonPushedFcn = createCallbackFcn(app, @UpdatelightingButtonPushed, true);
            app.UpdatelightingButton.Position = [426 19 100 23];
            app.UpdatelightingButton.Text = 'Update lighting';

            % Create TransparencySpinnerLabel
            app.TransparencySpinnerLabel = uilabel(app.PSFThreeDimUIFigure);
            app.TransparencySpinnerLabel.HorizontalAlignment = 'right';
            app.TransparencySpinnerLabel.Position = [548 19 81 22];
            app.TransparencySpinnerLabel.Text = 'Transparency';

            % Create TransparencySpinner
            app.TransparencySpinner = uispinner(app.PSFThreeDimUIFigure);
            app.TransparencySpinner.Step = 0.1;
            app.TransparencySpinner.Limits = [0.1 1];
            app.TransparencySpinner.ValueChangedFcn = createCallbackFcn(app, @TransparencySpinnerValueChanged, true);
            app.TransparencySpinner.BusyAction = 'cancel';
            app.TransparencySpinner.Position = [642 19 61 22];
            app.TransparencySpinner.Value = 0.2;

            % Show the figure after all components are created
            app.PSFThreeDimUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowPlotPSFThreeDim_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.PSFThreeDimUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.PSFThreeDimUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PSFThreeDimUIFigure)
        end
    end
end