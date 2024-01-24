classdef WindowGenerateDataset_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GeneratedatasetUIFigure       matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        RandomizedparametersLabel     matlab.ui.control.Label
        MaxLabel                      matlab.ui.control.Label
        MinLabel                      matlab.ui.control.Label
        RotationalfreedomCheckBox     matlab.ui.control.CheckBox
        DipoleorientationCheckBox     matlab.ui.control.CheckBox
        SaveimagesCheckBox            matlab.ui.control.CheckBox
        OutputFolder                  matlab.ui.control.Label
        OutputpathLabel               matlab.ui.control.Label
        SelectfolderButton            matlab.ui.control.Button
        NumberofgeneratedPSFsSpinner  matlab.ui.control.Spinner
        NumberofgeneratedPSFsSpinnerLabel  matlab.ui.control.Label
        SelectparametersandsetrangesLabel  matlab.ui.control.Label
        EmissionwavelengthCheckBox    matlab.ui.control.CheckBox
        WavelengthLowerLimitSpinner   matlab.ui.control.Spinner
        WavelengthUpperLimitSpinner   matlab.ui.control.Spinner
        PhotonnumberCheckBox          matlab.ui.control.CheckBox
        PhotonsLowerLimitSpinner      matlab.ui.control.Spinner
        PhotonsUpperLimitSpinner      matlab.ui.control.Spinner
        xpositionCheckBox             matlab.ui.control.CheckBox
        xPositionLowerLimitSpinner    matlab.ui.control.Spinner
        xPositionUpperLimitSpinner    matlab.ui.control.Spinner
        ypositionCheckBox             matlab.ui.control.CheckBox
        yPositionLowerLimitSpinner    matlab.ui.control.Spinner
        yPositionUpperLimitSpinner    matlab.ui.control.Spinner
        zpositionCheckBox             matlab.ui.control.CheckBox
        zPositionLowerLimitSpinner    matlab.ui.control.Spinner
        zPositionUpperLimitSpinner    matlab.ui.control.Spinner
        DefocusCheckBox               matlab.ui.control.CheckBox
        DefocusLowerLimitSpinner      matlab.ui.control.Spinner
        DefocusUpperLimitSpinner      matlab.ui.control.Spinner
        GeneraterandomizeddataButton  matlab.ui.control.Button
    end

    
    properties (Access = private)
        CallingApp
        outputFolderPath = ''
        outputFolderPathPsfDataSubfolder = ''
        outputFolderPathPsfImagesSubfolder = ''
    end
    
    methods (Access = private)
        
        function randomValues = getRandomParameter(app, n, lowerLimit, upperLimit, checkBoxValue, defaultValue)
            if checkBoxValue
                randomValues = lowerLimit + (upperLimit-lowerLimit) * rand(n,1);
            else
                % Repeat default parameter n times
                randomValues = repmat(defaultValue, n, 1);
            end
        end

        function randomValues = getRandomIntegerParameter(app, n, lowerLimit, upperLimit, checkBoxValue, defaultValue)
            if checkBoxValue
                randomValues = randi([lowerLimit, upperLimit], n, 1);
            else
                % Repeat default parameter n times
                randomValues = repmat(defaultValue, n, 1);
            end
        end

        function [inclination, azimuth] = getRandomDipole(app, n, checkBoxValue)
            if checkBoxValue
                inclination = asin(rand(n,1)); % angle between 0 and pi/2
                azimuth = 2*pi*rand(n,1); % angle between 0 and 2*pi
            else
                % Repeat default parameter n times
                inclination = repmat(app.CallingApp.theta.Value*pi/180, n, 1);
                azimuth = repmat(app.CallingApp.phi.Value*pi/180, n, 1);
            end
        end
    end

    methods (Access = private)
        function savePsfAsImage(app, psfValues, pngFullFilePath)
            f = figure('visible', 'off');
            imagesc(psfValues);
            axis equal; axis tight; axis off;
            colormap('viridis');
            exportgraphics(gca, pngFullFilePath, 'Resolution', 300);
            close(f);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.CallingApp = mainapp; % Store main app object

            if strcmp(app.CallingApp.DipolerotationButtonGroup.SelectedObject.Text, 'partially rotating')
                app.RotationalfreedomCheckBox.Visible = "on";
            else
                app.RotationalfreedomCheckBox.Visible = "off";
            end
        end

        % Button pushed function: GeneraterandomizeddataButton
        function GeneraterandomizeddataButtonPushed(app, event)
            if isempty(app.outputFolderPath)
                % Select output path before proceeding
                app.SelectfolderButtonPushed()
            end
            
            % Create output folders
            % Create subfolder with timestamp
            t = datetime('now','TimeZone','local','Format','yyyy-MM-dd_HH-mm-ss');
            timestamp = datestr(t, 'yyyy-MM-dd_HH-mm-ss');
            parentFolder = fullfile(app.outputFolderPath, ['psf_dataset_', timestamp]);
            mkdir(parentFolder);

            subfolderDataName = 'psf_data';
            subfolderDataPath = fullfile(parentFolder, subfolderDataName);
            if ~exist(subfolderDataPath, 'dir')
                mkdir(subfolderDataPath);
            end
            app.outputFolderPathPsfDataSubfolder = subfolderDataPath;
            
            if app.SaveimagesCheckBox.Value
                subfolderImagesName = 'psf_images';
                subfolderImagesPath = fullfile(parentFolder, subfolderImagesName);
                if ~exist(subfolderImagesPath, 'dir')
                    mkdir(subfolderImagesPath);
                end
                app.outputFolderPathPsfImagesSubfolder = subfolderImagesPath;
            end


            % Progress bar
            progress = uiprogressdlg(app.GeneratedatasetUIFigure, 'Title', 'Running...',...
                                     'Message', 'Creating dataset...', 'Cancelable', 'on');
            
            % Get all parameters
            par = app.CallingApp.getParameters();
            
            n = app.NumberofgeneratedPSFsSpinner.Value;
            parRandomWavelengths = app.getRandomParameter(n, app.WavelengthLowerLimitSpinner.Value, app.WavelengthUpperLimitSpinner.Value, app.EmissionwavelengthCheckBox.Value, app.CallingApp.EmissionWavelengthSpinner.Value);
            parRandomPhotonNumbers = app.getRandomIntegerParameter(n, app.PhotonsLowerLimitSpinner.Value, app.PhotonsUpperLimitSpinner.Value, app.PhotonnumberCheckBox.Value, app.CallingApp.PhotonnumberSpinner.Value);
            parRandomXPositions = app.getRandomParameter(n, app.xPositionLowerLimitSpinner.Value, app.xPositionUpperLimitSpinner.Value, app.xpositionCheckBox.Value, app.CallingApp.xpositionSpinner.Value);
            parRandomYPositions = app.getRandomParameter(n, app.yPositionLowerLimitSpinner.Value, app.yPositionUpperLimitSpinner.Value, app.ypositionCheckBox.Value, app.CallingApp.ypositionSpinner.Value);
            parRandomZPositions = app.getRandomParameter(n, app.zPositionLowerLimitSpinner.Value, app.zPositionUpperLimitSpinner.Value, app.zpositionCheckBox.Value, app.CallingApp.zpositionSpinner.Value);
            parRandomDefocus = app.getRandomParameter(n, app.DefocusLowerLimitSpinner.Value, app.DefocusUpperLimitSpinner.Value, app.DefocusCheckBox.Value, app.CallingApp.defocus.Value);
            [parRandomInclinationAngles, parRandomAzimuthAngles] = app.getRandomDipole(n, app.DipoleorientationCheckBox.Value);
            if strcmp(app.CallingApp.DipolerotationButtonGroup.SelectedObject.Text, 'partially rotating') && app.RotationalfreedomCheckBox.Value
                parRotationalConstraints = rand(n,1);
            else
                parRotationalConstraints = repmat(app.CallingApp.RotationalfreedomSpinner.Value, n, 1);
            end

            % Generate data
            numDigits = floor(log10(n)) + 1;
            for k = 1:n
                % Update progress bar
                progress.Value = k / n;
                progress.Message = sprintf('Creating data %d of %d', k, n);
                if progress.CancelRequested
                    break;
                end

                % Get random parameters
                par.wavelength = Length(parRandomWavelengths(k), 'nm');
                par.nPhotons = parRandomPhotonNumbers(k);
                par.position = Length([parRandomXPositions(k), parRandomYPositions(k), parRandomZPositions(k)], 'nm');
                par.defocus = Length([parRandomDefocus(k)], 'nm');
                par.dipole = Dipole(parRandomInclinationAngles(k),parRandomAzimuthAngles(k));
                if strcmp(app.CallingApp.DipolerotationButtonGroup.SelectedObject.Text, 'partially rotating') && app.RotationalfreedomCheckBox.Value
                    par.rotationalConstraint = parRotationalConstraints(k);
                end
                
                % Simulate PSF
                switch app.CallingApp.DipolerotationButtonGroup.SelectedObject.Text
                    case 'fixed'
                        psf = PSF(par);
                    case 'freely rotating'
                        psf = IsotropicPSF(par);
                    case 'partially rotating'
                        psf = PartiallyIsotropicPSF(par);
                    otherwise
                        error('Invalid input value for dipole rotation!')
                end
                psfValues = psf.image;

                % Save PSF
                try
                    % PSF data
                    filename = sprintf('psf_%0*d.mat', numDigits, k);
                    fullFilePath = fullfile(app.outputFolderPathPsfDataSubfolder, filename);
                    save(fullFilePath, 'psfValues');
                    
                    if app.SaveimagesCheckBox.Value
                        % PSF image
                        pngFilename = sprintf('psf_%0*d.png', numDigits, k);
                        pngFullFilePath = fullfile(app.outputFolderPathPsfImagesSubfolder, pngFilename);                
                        app.savePsfAsImage(psfValues, pngFullFilePath)
                    end
                catch ME
                    % Error handling
                    warndlg(['Error saving file: ', ME.message],'Warning');
                end
            end

            % Save fixed parameters
            fixedParameters = par;
            fieldsToRemove = {'wavelength', 'nPhotons', 'position', 'defocus', 'dipole', 'rotationalConstraint'};
            fixedParameters = rmfield(fixedParameters, fieldsToRemove);
            fileName = 'fixedParameters.mat';
            fullPath = fullfile(parentFolder, fileName);
            save(fullPath, 'fixedParameters');
            
            % Save table of changing parameters
            if strcmp(app.CallingApp.DipolerotationButtonGroup.SelectedObject.Text, 'partially rotating')
                randomizedParameterTable = table(parRandomWavelengths, parRandomPhotonNumbers, ...
                                             parRandomXPositions, parRandomYPositions, parRandomZPositions, parRandomDefocus, ...
                                             parRandomInclinationAngles*180/pi, parRandomAzimuthAngles*180/pi, parRotationalConstraints, ...
                                             'VariableNames', {'wavelength', 'nPhotons', 'x-position', 'y-position', 'z-position', 'defocus', ...
                                             'inclinationAngle', 'azimuthalAngle', 'rotationalConstraint'});
            else
                randomizedParameterTable = table(parRandomWavelengths, parRandomPhotonNumbers, ...
                                             parRandomXPositions, parRandomYPositions, parRandomZPositions, parRandomDefocus, ...
                                             parRandomInclinationAngles*180/pi, parRandomAzimuthAngles*180/pi, ...
                                             'VariableNames', {'wavelength', 'nPhotons', 'x-position', 'y-position', 'z-position', 'defocus', ...
                                             'inclinationAngle', 'azimuthalAngle'});
            end
            fileName = 'randomizedParameters.csv';
            fullFilePath = fullfile(parentFolder, fileName);
            writetable(randomizedParameterTable, fullFilePath);

            % Open output folder
            switch computer % Determine the operating system
                case 'PCWIN64', commandStr = ['explorer "', app.outputFolderPath, '"'];
                case 'GLNXA64', commandStr = ['xdg-open "', app.outputFolderPath, '"'];
                case 'MACI64', commandStr = ['open "', app.outputFolderPath, '"'];
                otherwise, error('Platform not supported')
            end
            system(commandStr);
        end

        % Button pushed function: SelectfolderButton
        function SelectfolderButtonPushed(app, event)
            folderpath = uigetdir;
            if folderpath == 0
                disp('User cancelled the folder selection.');
            else
                app.outputFolderPath = folderpath;
                app.OutputFolder.Text = folderpath;
            end
        end

        % Close request function: GeneratedatasetUIFigure
        function GeneratedatasetUIFigureCloseRequest(app, event)
            app.CallingApp.SwitchSubwindowGenerateDataset = 0;
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GeneratedatasetUIFigure and hide until all components are created
            app.GeneratedatasetUIFigure = uifigure('Visible', 'off');
            app.GeneratedatasetUIFigure.Position = [100 100 391 511];
            app.GeneratedatasetUIFigure.Name = 'Generate dataset';
            app.GeneratedatasetUIFigure.CloseRequestFcn = createCallbackFcn(app, @GeneratedatasetUIFigureCloseRequest, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.GeneratedatasetUIFigure);
            app.GridLayout.ColumnWidth = {'1x', 100, 100};
            app.GridLayout.RowHeight = {22, 15, 22, 22, 22, 22, 22, 22, 22, 22, 10, 22, 22, 22, 22, '1x'};

            % Create GeneraterandomizeddataButton
            app.GeneraterandomizeddataButton = uibutton(app.GridLayout, 'push');
            app.GeneraterandomizeddataButton.ButtonPushedFcn = createCallbackFcn(app, @GeneraterandomizeddataButtonPushed, true);
            app.GeneraterandomizeddataButton.Layout.Row = 16;
            app.GeneraterandomizeddataButton.Layout.Column = [1 3];
            app.GeneraterandomizeddataButton.Text = 'Generate randomized data';

            % Create DefocusUpperLimitSpinner
            app.DefocusUpperLimitSpinner = uispinner(app.GridLayout);
            app.DefocusUpperLimitSpinner.Step = 100;
            app.DefocusUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.DefocusUpperLimitSpinner.Layout.Row = 8;
            app.DefocusUpperLimitSpinner.Layout.Column = 3;
            app.DefocusUpperLimitSpinner.Value = 1000;

            % Create DefocusLowerLimitSpinner
            app.DefocusLowerLimitSpinner = uispinner(app.GridLayout);
            app.DefocusLowerLimitSpinner.Step = 100;
            app.DefocusLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.DefocusLowerLimitSpinner.Layout.Row = 8;
            app.DefocusLowerLimitSpinner.Layout.Column = 2;
            app.DefocusLowerLimitSpinner.Value = -1000;

            % Create DefocusCheckBox
            app.DefocusCheckBox = uicheckbox(app.GridLayout);
            app.DefocusCheckBox.Text = 'Defocus';
            app.DefocusCheckBox.Layout.Row = 8;
            app.DefocusCheckBox.Layout.Column = 1;
            app.DefocusCheckBox.Value = true;

            % Create zPositionUpperLimitSpinner
            app.zPositionUpperLimitSpinner = uispinner(app.GridLayout);
            app.zPositionUpperLimitSpinner.Step = 10;
            app.zPositionUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.zPositionUpperLimitSpinner.Layout.Row = 7;
            app.zPositionUpperLimitSpinner.Layout.Column = 3;
            app.zPositionUpperLimitSpinner.Value = 200;

            % Create zPositionLowerLimitSpinner
            app.zPositionLowerLimitSpinner = uispinner(app.GridLayout);
            app.zPositionLowerLimitSpinner.Step = 10;
            app.zPositionLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.zPositionLowerLimitSpinner.Layout.Row = 7;
            app.zPositionLowerLimitSpinner.Layout.Column = 2;
            app.zPositionLowerLimitSpinner.Value = -200;

            % Create zpositionCheckBox
            app.zpositionCheckBox = uicheckbox(app.GridLayout);
            app.zpositionCheckBox.Text = 'z-position';
            app.zpositionCheckBox.Layout.Row = 7;
            app.zpositionCheckBox.Layout.Column = 1;
            app.zpositionCheckBox.Value = true;

            % Create yPositionUpperLimitSpinner
            app.yPositionUpperLimitSpinner = uispinner(app.GridLayout);
            app.yPositionUpperLimitSpinner.Step = 10;
            app.yPositionUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.yPositionUpperLimitSpinner.Layout.Row = 6;
            app.yPositionUpperLimitSpinner.Layout.Column = 3;
            app.yPositionUpperLimitSpinner.Value = 200;

            % Create yPositionLowerLimitSpinner
            app.yPositionLowerLimitSpinner = uispinner(app.GridLayout);
            app.yPositionLowerLimitSpinner.Step = 10;
            app.yPositionLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.yPositionLowerLimitSpinner.Layout.Row = 6;
            app.yPositionLowerLimitSpinner.Layout.Column = 2;
            app.yPositionLowerLimitSpinner.Value = -200;

            % Create ypositionCheckBox
            app.ypositionCheckBox = uicheckbox(app.GridLayout);
            app.ypositionCheckBox.Text = 'y-position';
            app.ypositionCheckBox.Layout.Row = 6;
            app.ypositionCheckBox.Layout.Column = 1;
            app.ypositionCheckBox.Value = true;

            % Create xPositionUpperLimitSpinner
            app.xPositionUpperLimitSpinner = uispinner(app.GridLayout);
            app.xPositionUpperLimitSpinner.Step = 10;
            app.xPositionUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.xPositionUpperLimitSpinner.Layout.Row = 5;
            app.xPositionUpperLimitSpinner.Layout.Column = 3;
            app.xPositionUpperLimitSpinner.Value = 200;

            % Create xPositionLowerLimitSpinner
            app.xPositionLowerLimitSpinner = uispinner(app.GridLayout);
            app.xPositionLowerLimitSpinner.Step = 10;
            app.xPositionLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.xPositionLowerLimitSpinner.Layout.Row = 5;
            app.xPositionLowerLimitSpinner.Layout.Column = 2;
            app.xPositionLowerLimitSpinner.Value = -200;

            % Create xpositionCheckBox
            app.xpositionCheckBox = uicheckbox(app.GridLayout);
            app.xpositionCheckBox.Text = 'x-position';
            app.xpositionCheckBox.Layout.Row = 5;
            app.xpositionCheckBox.Layout.Column = 1;
            app.xpositionCheckBox.Value = true;

            % Create PhotonsUpperLimitSpinner
            app.PhotonsUpperLimitSpinner = uispinner(app.GridLayout);
            app.PhotonsUpperLimitSpinner.Step = 1000;
            app.PhotonsUpperLimitSpinner.Limits = [10 1000000];
            app.PhotonsUpperLimitSpinner.ValueDisplayFormat = '%d';
            app.PhotonsUpperLimitSpinner.Layout.Row = 4;
            app.PhotonsUpperLimitSpinner.Layout.Column = 3;
            app.PhotonsUpperLimitSpinner.Value = 100000;

            % Create PhotonsLowerLimitSpinner
            app.PhotonsLowerLimitSpinner = uispinner(app.GridLayout);
            app.PhotonsLowerLimitSpinner.Step = 100;
            app.PhotonsLowerLimitSpinner.Limits = [10 1000000];
            app.PhotonsLowerLimitSpinner.ValueDisplayFormat = '%d';
            app.PhotonsLowerLimitSpinner.Layout.Row = 4;
            app.PhotonsLowerLimitSpinner.Layout.Column = 2;
            app.PhotonsLowerLimitSpinner.Value = 100;

            % Create PhotonnumberCheckBox
            app.PhotonnumberCheckBox = uicheckbox(app.GridLayout);
            app.PhotonnumberCheckBox.Text = 'Photon number';
            app.PhotonnumberCheckBox.Layout.Row = 4;
            app.PhotonnumberCheckBox.Layout.Column = 1;
            app.PhotonnumberCheckBox.Value = true;

            % Create WavelengthUpperLimitSpinner
            app.WavelengthUpperLimitSpinner = uispinner(app.GridLayout);
            app.WavelengthUpperLimitSpinner.Step = 20;
            app.WavelengthUpperLimitSpinner.Limits = [200 1500];
            app.WavelengthUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.WavelengthUpperLimitSpinner.Layout.Row = 3;
            app.WavelengthUpperLimitSpinner.Layout.Column = 3;
            app.WavelengthUpperLimitSpinner.Value = 700;

            % Create WavelengthLowerLimitSpinner
            app.WavelengthLowerLimitSpinner = uispinner(app.GridLayout);
            app.WavelengthLowerLimitSpinner.Step = 20;
            app.WavelengthLowerLimitSpinner.Limits = [200 1500];
            app.WavelengthLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.WavelengthLowerLimitSpinner.Layout.Row = 3;
            app.WavelengthLowerLimitSpinner.Layout.Column = 2;
            app.WavelengthLowerLimitSpinner.Value = 450;

            % Create EmissionwavelengthCheckBox
            app.EmissionwavelengthCheckBox = uicheckbox(app.GridLayout);
            app.EmissionwavelengthCheckBox.Text = 'Emission wavelength';
            app.EmissionwavelengthCheckBox.Layout.Row = 3;
            app.EmissionwavelengthCheckBox.Layout.Column = 1;

            % Create SelectparametersandsetrangesLabel
            app.SelectparametersandsetrangesLabel = uilabel(app.GridLayout);
            app.SelectparametersandsetrangesLabel.FontSize = 13;
            app.SelectparametersandsetrangesLabel.FontWeight = 'bold';
            app.SelectparametersandsetrangesLabel.Layout.Row = 1;
            app.SelectparametersandsetrangesLabel.Layout.Column = [1 3];
            app.SelectparametersandsetrangesLabel.Text = 'Select parameters and set ranges';

            % Create NumberofgeneratedPSFsSpinnerLabel
            app.NumberofgeneratedPSFsSpinnerLabel = uilabel(app.GridLayout);
            app.NumberofgeneratedPSFsSpinnerLabel.Layout.Row = 12;
            app.NumberofgeneratedPSFsSpinnerLabel.Layout.Column = 1;
            app.NumberofgeneratedPSFsSpinnerLabel.Text = 'Number of generated PSFs';

            % Create NumberofgeneratedPSFsSpinner
            app.NumberofgeneratedPSFsSpinner = uispinner(app.GridLayout);
            app.NumberofgeneratedPSFsSpinner.Step = 100;
            app.NumberofgeneratedPSFsSpinner.Limits = [1 10000];
            app.NumberofgeneratedPSFsSpinner.RoundFractionalValues = 'on';
            app.NumberofgeneratedPSFsSpinner.ValueDisplayFormat = '%d';
            app.NumberofgeneratedPSFsSpinner.Layout.Row = 12;
            app.NumberofgeneratedPSFsSpinner.Layout.Column = 2;
            app.NumberofgeneratedPSFsSpinner.Value = 100;

            % Create SelectfolderButton
            app.SelectfolderButton = uibutton(app.GridLayout, 'push');
            app.SelectfolderButton.ButtonPushedFcn = createCallbackFcn(app, @SelectfolderButtonPushed, true);
            app.SelectfolderButton.Layout.Row = 13;
            app.SelectfolderButton.Layout.Column = 3;
            app.SelectfolderButton.Text = 'Select folder';

            % Create OutputpathLabel
            app.OutputpathLabel = uilabel(app.GridLayout);
            app.OutputpathLabel.Layout.Row = 13;
            app.OutputpathLabel.Layout.Column = [1 2];
            app.OutputpathLabel.Text = 'Output path:';

            % Create OutputFolder
            app.OutputFolder = uilabel(app.GridLayout);
            app.OutputFolder.Layout.Row = 14;
            app.OutputFolder.Layout.Column = [1 3];
            app.OutputFolder.Text = '';

            % Create SaveimagesCheckBox
            app.SaveimagesCheckBox = uicheckbox(app.GridLayout);
            app.SaveimagesCheckBox.Text = 'Save images';
            app.SaveimagesCheckBox.Layout.Row = 15;
            app.SaveimagesCheckBox.Layout.Column = 1;

            % Create DipoleorientationCheckBox
            app.DipoleorientationCheckBox = uicheckbox(app.GridLayout);
            app.DipoleorientationCheckBox.Tooltip = {'Random orientational distribution on sphere'};
            app.DipoleorientationCheckBox.Text = 'Dipole orientation';
            app.DipoleorientationCheckBox.Layout.Row = 9;
            app.DipoleorientationCheckBox.Layout.Column = 1;
            app.DipoleorientationCheckBox.Value = true;

            % Create RotationalfreedomCheckBox
            app.RotationalfreedomCheckBox = uicheckbox(app.GridLayout);
            app.RotationalfreedomCheckBox.Tooltip = {'Range 0-1'};
            app.RotationalfreedomCheckBox.Text = 'Rotational freedom';
            app.RotationalfreedomCheckBox.Layout.Row = 10;
            app.RotationalfreedomCheckBox.Layout.Column = 1;
            app.RotationalfreedomCheckBox.Value = true;

            % Create MinLabel
            app.MinLabel = uilabel(app.GridLayout);
            app.MinLabel.HorizontalAlignment = 'center';
            app.MinLabel.VerticalAlignment = 'bottom';
            app.MinLabel.Layout.Row = 2;
            app.MinLabel.Layout.Column = 2;
            app.MinLabel.Text = 'Min.';

            % Create MaxLabel
            app.MaxLabel = uilabel(app.GridLayout);
            app.MaxLabel.HorizontalAlignment = 'center';
            app.MaxLabel.VerticalAlignment = 'bottom';
            app.MaxLabel.Layout.Row = 2;
            app.MaxLabel.Layout.Column = 3;
            app.MaxLabel.Text = 'Max.';

            % Create RandomizedparametersLabel
            app.RandomizedparametersLabel = uilabel(app.GridLayout);
            app.RandomizedparametersLabel.FontWeight = 'bold';
            app.RandomizedparametersLabel.Tooltip = {'If parameter is not selected, its value will be taken from the main app'};
            app.RandomizedparametersLabel.Layout.Row = 2;
            app.RandomizedparametersLabel.Layout.Column = 1;
            app.RandomizedparametersLabel.Text = 'Randomized parameters';

            % Show the figure after all components are created
            app.GeneratedatasetUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowGenerateDataset_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.GeneratedatasetUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.GeneratedatasetUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GeneratedatasetUIFigure)
        end
    end
end