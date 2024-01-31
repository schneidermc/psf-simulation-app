classdef WindowPolarizedEmission_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PolarizedEmissionUIFigure  matlab.ui.Figure
        Toolbar                    matlab.ui.container.Toolbar
        SaveXPolarized             matlab.ui.container.toolbar.PushTool
        SaveYPolarized             matlab.ui.container.toolbar.PushTool
        ShareColorbar              matlab.ui.container.toolbar.PushTool
        UIAxesPolarizedEmission_y  matlab.ui.control.UIAxes
        UIAxesPolarizedEmission_x  matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        CallingApp   % Main app object
    end


    
    methods (Access = public)
        function initializePlot(app)
            imagesc(app.UIAxesPolarizedEmission_x, zeros(app.CallingApp.nPixels))
            axis(app.UIAxesPolarizedEmission_x, 'equal');
            axis(app.UIAxesPolarizedEmission_x, 'tight');
            title(app.UIAxesPolarizedEmission_x, 'x-polarized')
            set(app.UIAxesPolarizedEmission_x,'visible','on')
            imagesc(app.UIAxesPolarizedEmission_y, zeros(app.CallingApp.nPixels))
            axis(app.UIAxesPolarizedEmission_y, 'equal');
            axis(app.UIAxesPolarizedEmission_y, 'tight');  
            title(app.UIAxesPolarizedEmission_y, 'y-polarized')
            set(app.UIAxesPolarizedEmission_y,'visible','on')
            colormap(app.PolarizedEmissionUIFigure, app.CallingApp.ColormapDropDown.Value)
        end

        function updatePlot(app, psf)
            if app.CallingApp.ShowPsf3DCheckBox.Value
                nSteps = app.CallingApp.Numberzsteps3DPSFEditField.Value - 1;
                midSlice = floor(nSteps/2);
                Ix = psf.Ix(:,:,midSlice);
                Iy = psf.Iy(:,:,midSlice);
                imagesc(app.UIAxesPolarizedEmission_x, Ix)
                imagesc(app.UIAxesPolarizedEmission_y, Iy)
            else
                imagesc(app.UIAxesPolarizedEmission_x, psf.Ix)
                imagesc(app.UIAxesPolarizedEmission_y, psf.Iy)
            end
            colorbar(app.UIAxesPolarizedEmission_y)
            switch app.ShareColorbar.Icon
                case 'brushActive.png'
                    limits = ([0 max(max(psf.Ix(:), psf.Iy(:)))]);
                    caxis(app.UIAxesPolarizedEmission_x, limits)
                    caxis(app.UIAxesPolarizedEmission_y, limits)
                    %colorbar(app.UIAxesPolarizedEmission_x,'off')
                case 'brushInactive.png'
                    colorbar(app.UIAxesPolarizedEmission_x)
                    caxis(app.UIAxesPolarizedEmission_x, 'auto')
                    caxis(app.UIAxesPolarizedEmission_y, 'auto')
            end
        end

        function setColormap(app,map)
            colormap(app.PolarizedEmissionUIFigure, map)
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.CallingApp = mainapp;
        end

        % Close request function: PolarizedEmissionUIFigure
        function PolarizedEmissionUIFigureCloseRequest(app, event)
            app.CallingApp.PolarizedemissionchannelsCheckBox.Value = 0;
            delete(app)
        end

        % Callback function
        function SavexpolarizedasMenuSelected(app, event)
            psf_x = getimage(app.UIAxesPolarizedEmission_x);
            startingFolder = userpath;
            defaultFileName = fullfile(startingFolder, 'psf_x.csv');
            [filename, folder] = uiputfile(defaultFileName, 'Specify a file');
            if filename == 0
              % User clicked the Cancel button.
              return;
            end
            % save file 
            filename = fullfile(folder,filename); 
            writematrix(psf_x, filename);
        end

        % Callback function
        function SaveypolarizedasMenuSelected(app, event)
            psf_y = getimage(app.UIAxesPolarizedEmission_y);
            startingFolder = userpath;
            defaultFileName = fullfile(startingFolder, 'psf_y.csv');
            [filename, folder] = uiputfile(defaultFileName, 'Specify a file');
            if filename == 0
              % User clicked the Cancel button.
              return;
            end
            % save file 
            filename = fullfile(folder,filename); 
            writematrix(psf_y, filename);
        end

        % Clicked callback: ShareColorbar
        function ShareColorbarClicked(app, event)
            switch app.ShareColorbar.Icon
                case 'brushInactive.png'
                    app.ShareColorbar.Icon = 'brushActive.png';
                case 'brushActive.png'
                    app.ShareColorbar.Icon = 'brushInactive.png';
            end
            switch app.CallingApp.SwitchMultipleFluorophores
                case 0 % 'Single'
                    psf = app.CallingApp.simulateAndDisplayPSF;
                    app.CallingApp.PlotPolarizedEmission.updatePlot(psf);
                case 1 % 'Multiple'
                    app.CallingApp.Fluorophores.updatePolarizedEmissionChannels();
            end
        end

        % Clicked callback: SaveXPolarized
        function SaveXPolarizedClicked(app, event)
            psf_x = getimage(app.UIAxesPolarizedEmission_x);
            startingFolder = userpath;
            filter = {'*.dat'; '*.xls'; '*.xlsx'; '*.csv'; '*.txt'; '*.mat'; '*.png'; '*.jpg'; '*.tif'};
            defaultFileName = fullfile(startingFolder, filter);
            [filename, folder] = uiputfile(defaultFileName, 'Specify a file', 'psf_x');
            if filename == 0
              % User clicked the Cancel button.
              return;
            end
            [~,~,ext] = fileparts(filename); 
            filename = fullfile(folder,filename);  
            switch ext 
                case {'.dat', '.xls', '.xlsx', '.csv', '.txt'}
                    writematrix(psf_x, filename);
                case '.mat'
                    save(filename,'psf_x');
                case {'.png', '.tif'}
                    imwrite( ind2rgb(im2uint8(mat2gray(psf_x)), colormap(app.PolarizedEmissionUIFigure, app.CallingApp.ColormapDropDown.Value)), filename)
                case '.jpg'
                    Nx = app.CallingApp.PixelsperlateralaxisEditField.Value; 
                    psf_x = interp2(1:Nx, (1:Nx)', psf_x, 1:0.02:Nx, (1:0.02:Nx)', 'nearest');
                    imwrite( ind2rgb(im2uint8(mat2gray(psf_x)), colormap(app.PolarizedEmissionUIFigure, app.CallingApp.ColormapDropDown.Value)), filename)
            end
        end

        % Clicked callback: SaveYPolarized
        function SaveYPolarizedClicked(app, event)
            psf_y = getimage(app.UIAxesPolarizedEmission_y);
            startingFolder = userpath;
            filter = {'*.dat'; '*.xls'; '*.xlsx'; '*.csv'; '*.txt'; '*.mat'; '*.png'; '*.jpg'; '*.tif'};
            defaultFileName = fullfile(startingFolder, filter);
            [filename, folder] = uiputfile(defaultFileName, 'Specify a file', 'psf_y');
            if filename == 0
              % User clicked the Cancel button.
              return;
            end
            [~,~,ext] = fileparts(filename); 
            filename = fullfile(folder,filename);  
            switch ext 
                case {'.dat', '.xls', '.xlsx', '.csv', '.txt'}
                    writematrix(psf_y, filename);
                case '.mat'
                    save(filename,'psf_y');
                case '.tif'
                    imwrite( ind2rgb(im2uint8(mat2gray(psf_y)), colormap(app.PolarizedEmissionUIFigure, app.CallingApp.ColormapDropDown.Value)), filename)
                case {'.png', '.jpg'}
                    Nx = app.CallingApp.PixelsperlateralaxisEditField.Value; 
                    psf_y = interp2(1:Nx, (1:Nx)', psf_y, 1:0.02:Nx, (1:0.02:Nx)', 'nearest');
                    imwrite( ind2rgb(im2uint8(mat2gray(psf_y)), colormap(app.PolarizedEmissionUIFigure, app.CallingApp.ColormapDropDown.Value)), filename)
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PolarizedEmissionUIFigure and hide until all components are created
            app.PolarizedEmissionUIFigure = uifigure('Visible', 'off');
            app.PolarizedEmissionUIFigure.Position = [100 100 601 293];
            app.PolarizedEmissionUIFigure.Name = 'Polarized emission';
            app.PolarizedEmissionUIFigure.CloseRequestFcn = createCallbackFcn(app, @PolarizedEmissionUIFigureCloseRequest, true);
            app.PolarizedEmissionUIFigure.BusyAction = 'cancel';
            app.PolarizedEmissionUIFigure.HandleVisibility = 'callback';

            % Create Toolbar
            app.Toolbar = uitoolbar(app.PolarizedEmissionUIFigure);

            % Create SaveXPolarized
            app.SaveXPolarized = uipushtool(app.Toolbar);
            app.SaveXPolarized.Tooltip = {'Save x-polarized'};
            app.SaveXPolarized.ClickedCallback = createCallbackFcn(app, @SaveXPolarizedClicked, true);
            app.SaveXPolarized.Icon = 'save.png';

            % Create SaveYPolarized
            app.SaveYPolarized = uipushtool(app.Toolbar);
            app.SaveYPolarized.Tooltip = {'Save y-polarized'};
            app.SaveYPolarized.ClickedCallback = createCallbackFcn(app, @SaveYPolarizedClicked, true);
            app.SaveYPolarized.Icon = 'save.png';

            % Create ShareColorbar
            app.ShareColorbar = uipushtool(app.Toolbar);
            app.ShareColorbar.Tooltip = {'Share Colorbar'};
            app.ShareColorbar.ClickedCallback = createCallbackFcn(app, @ShareColorbarClicked, true);
            app.ShareColorbar.Icon = 'brushInactive.png';

            % Create UIAxesPolarizedEmission_x
            app.UIAxesPolarizedEmission_x = uiaxes(app.PolarizedEmissionUIFigure);
            zlabel(app.UIAxesPolarizedEmission_x, 'Z')
            app.UIAxesPolarizedEmission_x.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesPolarizedEmission_x.XTick = [];
            app.UIAxesPolarizedEmission_x.YTick = [];
            app.UIAxesPolarizedEmission_x.FontSize = 15;
            app.UIAxesPolarizedEmission_x.Position = [19 12 274 263];

            % Create UIAxesPolarizedEmission_y
            app.UIAxesPolarizedEmission_y = uiaxes(app.PolarizedEmissionUIFigure);
            zlabel(app.UIAxesPolarizedEmission_y, 'Z')
            app.UIAxesPolarizedEmission_y.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesPolarizedEmission_y.XTick = [];
            app.UIAxesPolarizedEmission_y.YTick = [];
            app.UIAxesPolarizedEmission_y.FontSize = 15;
            app.UIAxesPolarizedEmission_y.Position = [299 12 274 263];

            % Show the figure after all components are created
            app.PolarizedEmissionUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowPolarizedEmission_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PolarizedEmissionUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PolarizedEmissionUIFigure)
        end
    end
end