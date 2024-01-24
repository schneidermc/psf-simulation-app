classdef WindowPolarizedEmission_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PolarizedEmissionUIFigure  matlab.ui.Figure
        OptionsMenu                matlab.ui.container.Menu
        ShareColorbarMenu          matlab.ui.container.Menu
        SavexpolarizedasMenu       matlab.ui.container.Menu
        SaveypolarizedasMenu       matlab.ui.container.Menu
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
                nSteps = app.CallingApp.NumberstepsEditField.Value - 1;
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
            if strcmp(app.ShareColorbarMenu.Checked,'on')
                limits = ([0 max(max(psf.Ix(:), psf.Iy(:)))]);
                caxis(app.UIAxesPolarizedEmission_x, limits)
                caxis(app.UIAxesPolarizedEmission_y, limits)
                %colorbar(app.UIAxesPolarizedEmission_x,'off')
            else
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

        % Menu selected function: ShareColorbarMenu
        function ShareColorbarMenuSelected(app, event)
            if strcmp(app.ShareColorbarMenu.Checked,'on')
                app.ShareColorbarMenu.Checked = 'off';
            else
                app.ShareColorbarMenu.Checked = 'on';
            end
            switch app.CallingApp.SwitchMultipleFluorophores
                case 0 % 'Single'
                    psf = app.CallingApp.simulateAndDisplayPSF;
                    app.CallingApp.PlotPolarizedEmission.updatePlot(psf);
                case 1 % 'Multiple'
                    app.CallingApp.Fluorophores.updatePolarizedEmissionChannels();
            end
        end

        % Menu selected function: SavexpolarizedasMenu
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

        % Menu selected function: SaveypolarizedasMenu
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

            % Create OptionsMenu
            app.OptionsMenu = uimenu(app.PolarizedEmissionUIFigure);
            app.OptionsMenu.Text = 'Options';

            % Create ShareColorbarMenu
            app.ShareColorbarMenu = uimenu(app.OptionsMenu);
            app.ShareColorbarMenu.MenuSelectedFcn = createCallbackFcn(app, @ShareColorbarMenuSelected, true);
            app.ShareColorbarMenu.Text = 'Share Colorbar';

            % Create SavexpolarizedasMenu
            app.SavexpolarizedasMenu = uimenu(app.OptionsMenu);
            app.SavexpolarizedasMenu.MenuSelectedFcn = createCallbackFcn(app, @SavexpolarizedasMenuSelected, true);
            app.SavexpolarizedasMenu.Text = 'Save x-polarized as';

            % Create SaveypolarizedasMenu
            app.SaveypolarizedasMenu = uimenu(app.OptionsMenu);
            app.SaveypolarizedasMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveypolarizedasMenuSelected, true);
            app.SaveypolarizedasMenu.Text = 'Save y-polarized as';

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