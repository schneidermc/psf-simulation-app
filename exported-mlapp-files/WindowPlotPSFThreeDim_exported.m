classdef WindowPlotPSFThreeDim_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PSFThreeDimUIFigure       matlab.ui.Figure
        Toolbar                   matlab.ui.container.Toolbar
        saveProjectionPlot        matlab.ui.container.toolbar.PushTool
        updateLighting            matlab.ui.container.toolbar.PushTool
        TransparencySpinner       matlab.ui.control.Spinner
        TransparencySpinnerLabel  matlab.ui.control.Label
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
            set(app.UIAxesPSF,'DataAspectRatio',[1 1 1])
            
            nPixels = app.CallingApp.nPixels;
            pixelSize = app.CallingApp.PixelsizeObjectSpaceSpinner.Value;
            xAxis = linspace(0, pixelSize*nPixels, nPixels) - pixelSize*nPixels/2;
            zAxis = app.CallingApp.parameters.defocus.inNanometer();
            
            % Show xz projection
            imagesc(app.UIAxesPSF, xAxis, zAxis, projection');
            zLimit = 1/2*app.CallingApp.zstepsize3DPSFEditField.Value*(app.CallingApp.Numberzsteps3DPSFEditField.Value-1);
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
            app.CallingApp.zstepsize3DPSFEditField.Visible = "off";
            app.CallingApp.zstepsize3DPSFEditField.Enable = "off";
            app.CallingApp.zstepsize3DPSFEditFieldLabel.Visible = "off";
            app.CallingApp.zstepsize3DPSFEditFieldLabel.Enable = "off";
            app.CallingApp.Numberzsteps3DPSFEditField.Visible = "off";
            app.CallingApp.Numberzsteps3DPSFEditField.Enable = "off";
            app.CallingApp.Numberzsteps3DPSFEditFieldLabel.Visible = "off";
            app.CallingApp.Numberzsteps3DPSFEditFieldLabel.Enable = "off";
            if app.CallingApp.ShowPsf2DCheckBox.Value == false
                app.CallingApp.ShowPSFCheckBox.Value = false;
                app.CallingApp.ShowPsf2DCheckBox.Enable = "off";
                app.CallingApp.ShowPsf3DCheckBox.Enable = "off";
            end
            delete(app)
        end

        % Clicked callback: saveProjectionPlot
        function saveProjectionPlotClicked(app, event)
            psfXZ = getimage(app.UIAxesPSF);
            startingFolder = userpath;
            filter = {'*.dat'; '*.xls'; '*.xlsx'; '*.csv'; '*.txt'; '*.mat'; '*.png'; '*.jpg'; '*.tif'};
            defaultFileName = fullfile(startingFolder, filter);
            [filename, folder] = uiputfile(defaultFileName, 'Specify a file', 'psfXZ');
            if filename == 0
              % User clicked the Cancel button.
              return;
            end
            [~,~,ext] = fileparts(filename); 
            filename = fullfile(folder,filename);  
            switch ext 
                case {'.dat', '.xls', '.xlsx', '.csv', '.txt'}
                    writematrix(psfXZ, filename);
                case '.mat'
                    save(filename,'psfXZ');
                case {'.tif'}
                    imwrite( ind2rgb(im2uint8(mat2gray(psfXZ)), colormap(app.PSFThreeDimUIFigure, app.CallingApp.ColormapDropDown.Value)), filename)
                case {'.png', '.jpg'}
                    Nx = app.CallingApp.PixelsperlateralaxisEditField.Value; 
                    psfXZ = interp2(1:Nx, (1:Nx)', psfXZ, 1:0.02:Nx, (1:0.02:Nx)', 'nearest');
                    imwrite( ind2rgb(im2uint8(mat2gray(psfXZ)), colormap(app.PSFThreeDimUIFigure, app.CallingApp.ColormapDropDown.Value)), filename)
            end
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

        % Value changed function: TransparencySpinner
        function TransparencySpinnerValueChanged(app, event)
            app.TransparencySpinner.Value = event.Value;
            app.isosurfaceTransparency = app.TransparencySpinner.Value;
            alpha(app.UIAxesIsoSurface, app.isosurfaceTransparency);
        end

        % Clicked callback: updateLighting
        function updateLightingClicked(app, event)
            camlight(app.cameraLight, 'right');
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

            % Create saveProjectionPlot
            app.saveProjectionPlot = uipushtool(app.Toolbar);
            app.saveProjectionPlot.Tooltip = {'Save xz-Projection'};
            app.saveProjectionPlot.ClickedCallback = createCallbackFcn(app, @saveProjectionPlotClicked, true);
            app.saveProjectionPlot.Icon = 'save.png';

            % Create updateLighting
            app.updateLighting = uipushtool(app.Toolbar);
            app.updateLighting.Tooltip = {'Update isosurface lighting'};
            app.updateLighting.ClickedCallback = createCallbackFcn(app, @updateLightingClicked, true);
            app.updateLighting.Icon = 'lightbulb.png';

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