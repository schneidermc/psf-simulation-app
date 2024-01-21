classdef WindowGenerateDataset_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GeneratedatasetUIFigure       matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        SaveimagesCheckBox            matlab.ui.control.CheckBox
        OutputFolder                  matlab.ui.control.Label
        OutputfolderLabel             matlab.ui.control.Label
        SelectfolderButton            matlab.ui.control.Button
        NumberofgeneratedPSFsSpinner  matlab.ui.control.Spinner
        NumberofgeneratedPSFsSpinnerLabel  matlab.ui.control.Label
        SetparameterrangesLabel       matlab.ui.control.Label
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
        end

        % Button pushed function: GeneraterandomizeddataButton
        function GeneraterandomizeddataButtonPushed(app, event)
            if isempty(app.outputFolderPath)
                % Select output folder before proceeding
                app.SelectfolderButtonPushed()
            end
            if app.SaveimagesCheckBox.Value && isempty(app.outputFolderPathPsfImagesSubfolder)
                subfolderImagesName = 'psf_images';
                subfolderImagesPath = fullfile(app.outputFolderPath, subfolderImagesName);
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
            fieldsToRemove = {'wavelength', 'nPhotons', 'position', 'defocus'};
            fixedParameters = rmfield(fixedParameters, fieldsToRemove);
            fileName = 'fixedParameters.mat';
            fullPath = fullfile(app.outputFolderPath, fileName);
            save(fullPath, 'fixedParameters');
            
            % Save table of changing parameters
            randomizedParameterTable = table(parRandomWavelengths, parRandomPhotonNumbers, parRandomXPositions, parRandomYPositions, parRandomZPositions, parRandomDefocus, ...
                                   'VariableNames', {'wavelength', 'nPhotons', 'x-position', 'y-position', 'z-position', 'defocus'});
            fileName = 'randomizedParameters.csv';
            fullFilePath = fullfile(app.outputFolderPath, fileName);
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

                subfolderDataName = 'psf_data';
                subfolderDataPath = fullfile(folderpath, subfolderDataName);
                if ~exist(subfolderDataPath, 'dir')
                    mkdir(subfolderDataPath);
                end
                app.outputFolderPathPsfDataSubfolder = subfolderDataPath;
                
                if app.SaveimagesCheckBox.Value
                    subfolderImagesName = 'psf_images';
                    subfolderImagesPath = fullfile(folderpath, subfolderImagesName);
                    if ~exist(subfolderImagesPath, 'dir')
                        mkdir(subfolderImagesPath);
                    end
                    app.outputFolderPathPsfImagesSubfolder = subfolderImagesPath;
                end
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GeneratedatasetUIFigure and hide until all components are created
            app.GeneratedatasetUIFigure = uifigure('Visible', 'off');
            app.GeneratedatasetUIFigure.Position = [100 100 391 400];
            app.GeneratedatasetUIFigure.Name = 'Generate dataset';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.GeneratedatasetUIFigure);
            app.GridLayout.ColumnWidth = {'1x', 100, 100};
            app.GridLayout.RowHeight = {22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, '1x'};

            % Create GeneraterandomizeddataButton
            app.GeneraterandomizeddataButton = uibutton(app.GridLayout, 'push');
            app.GeneraterandomizeddataButton.ButtonPushedFcn = createCallbackFcn(app, @GeneraterandomizeddataButtonPushed, true);
            app.GeneraterandomizeddataButton.Layout.Row = 12;
            app.GeneraterandomizeddataButton.Layout.Column = [1 3];
            app.GeneraterandomizeddataButton.Text = 'Generate randomized data';

            % Create DefocusUpperLimitSpinner
            app.DefocusUpperLimitSpinner = uispinner(app.GridLayout);
            app.DefocusUpperLimitSpinner.Step = 100;
            app.DefocusUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.DefocusUpperLimitSpinner.Layout.Row = 7;
            app.DefocusUpperLimitSpinner.Layout.Column = 3;
            app.DefocusUpperLimitSpinner.Value = 1000;

            % Create DefocusLowerLimitSpinner
            app.DefocusLowerLimitSpinner = uispinner(app.GridLayout);
            app.DefocusLowerLimitSpinner.Step = 100;
            app.DefocusLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.DefocusLowerLimitSpinner.Layout.Row = 7;
            app.DefocusLowerLimitSpinner.Layout.Column = 2;
            app.DefocusLowerLimitSpinner.Value = -1000;

            % Create DefocusCheckBox
            app.DefocusCheckBox = uicheckbox(app.GridLayout);
            app.DefocusCheckBox.Text = 'Defocus';
            app.DefocusCheckBox.Layout.Row = 7;
            app.DefocusCheckBox.Layout.Column = 1;

            % Create zPositionUpperLimitSpinner
            app.zPositionUpperLimitSpinner = uispinner(app.GridLayout);
            app.zPositionUpperLimitSpinner.Step = 10;
            app.zPositionUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.zPositionUpperLimitSpinner.Layout.Row = 6;
            app.zPositionUpperLimitSpinner.Layout.Column = 3;
            app.zPositionUpperLimitSpinner.Value = 200;

            % Create zPositionLowerLimitSpinner
            app.zPositionLowerLimitSpinner = uispinner(app.GridLayout);
            app.zPositionLowerLimitSpinner.Step = 10;
            app.zPositionLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.zPositionLowerLimitSpinner.Layout.Row = 6;
            app.zPositionLowerLimitSpinner.Layout.Column = 2;
            app.zPositionLowerLimitSpinner.Value = -200;

            % Create zpositionCheckBox
            app.zpositionCheckBox = uicheckbox(app.GridLayout);
            app.zpositionCheckBox.Text = 'z-position';
            app.zpositionCheckBox.Layout.Row = 6;
            app.zpositionCheckBox.Layout.Column = 1;

            % Create yPositionUpperLimitSpinner
            app.yPositionUpperLimitSpinner = uispinner(app.GridLayout);
            app.yPositionUpperLimitSpinner.Step = 10;
            app.yPositionUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.yPositionUpperLimitSpinner.Layout.Row = 5;
            app.yPositionUpperLimitSpinner.Layout.Column = 3;
            app.yPositionUpperLimitSpinner.Value = 200;

            % Create yPositionLowerLimitSpinner
            app.yPositionLowerLimitSpinner = uispinner(app.GridLayout);
            app.yPositionLowerLimitSpinner.Step = 10;
            app.yPositionLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.yPositionLowerLimitSpinner.Layout.Row = 5;
            app.yPositionLowerLimitSpinner.Layout.Column = 2;
            app.yPositionLowerLimitSpinner.Value = -200;

            % Create ypositionCheckBox
            app.ypositionCheckBox = uicheckbox(app.GridLayout);
            app.ypositionCheckBox.Text = 'y-position';
            app.ypositionCheckBox.Layout.Row = 5;
            app.ypositionCheckBox.Layout.Column = 1;

            % Create xPositionUpperLimitSpinner
            app.xPositionUpperLimitSpinner = uispinner(app.GridLayout);
            app.xPositionUpperLimitSpinner.Step = 10;
            app.xPositionUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.xPositionUpperLimitSpinner.Layout.Row = 4;
            app.xPositionUpperLimitSpinner.Layout.Column = 3;
            app.xPositionUpperLimitSpinner.Value = 200;

            % Create xPositionLowerLimitSpinner
            app.xPositionLowerLimitSpinner = uispinner(app.GridLayout);
            app.xPositionLowerLimitSpinner.Step = 10;
            app.xPositionLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.xPositionLowerLimitSpinner.Layout.Row = 4;
            app.xPositionLowerLimitSpinner.Layout.Column = 2;
            app.xPositionLowerLimitSpinner.Value = -200;

            % Create xpositionCheckBox
            app.xpositionCheckBox = uicheckbox(app.GridLayout);
            app.xpositionCheckBox.Text = 'x-position';
            app.xpositionCheckBox.Layout.Row = 4;
            app.xpositionCheckBox.Layout.Column = 1;

            % Create PhotonsUpperLimitSpinner
            app.PhotonsUpperLimitSpinner = uispinner(app.GridLayout);
            app.PhotonsUpperLimitSpinner.Step = 1000;
            app.PhotonsUpperLimitSpinner.Limits = [10 1000000];
            app.PhotonsUpperLimitSpinner.ValueDisplayFormat = '%d';
            app.PhotonsUpperLimitSpinner.Layout.Row = 3;
            app.PhotonsUpperLimitSpinner.Layout.Column = 3;
            app.PhotonsUpperLimitSpinner.Value = 100000;

            % Create PhotonsLowerLimitSpinner
            app.PhotonsLowerLimitSpinner = uispinner(app.GridLayout);
            app.PhotonsLowerLimitSpinner.Step = 100;
            app.PhotonsLowerLimitSpinner.Limits = [10 1000000];
            app.PhotonsLowerLimitSpinner.ValueDisplayFormat = '%d';
            app.PhotonsLowerLimitSpinner.Layout.Row = 3;
            app.PhotonsLowerLimitSpinner.Layout.Column = 2;
            app.PhotonsLowerLimitSpinner.Value = 100;

            % Create PhotonnumberCheckBox
            app.PhotonnumberCheckBox = uicheckbox(app.GridLayout);
            app.PhotonnumberCheckBox.Text = 'Photon number';
            app.PhotonnumberCheckBox.Layout.Row = 3;
            app.PhotonnumberCheckBox.Layout.Column = 1;

            % Create WavelengthUpperLimitSpinner
            app.WavelengthUpperLimitSpinner = uispinner(app.GridLayout);
            app.WavelengthUpperLimitSpinner.Step = 20;
            app.WavelengthUpperLimitSpinner.Limits = [200 1500];
            app.WavelengthUpperLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.WavelengthUpperLimitSpinner.Layout.Row = 2;
            app.WavelengthUpperLimitSpinner.Layout.Column = 3;
            app.WavelengthUpperLimitSpinner.Value = 700;

            % Create WavelengthLowerLimitSpinner
            app.WavelengthLowerLimitSpinner = uispinner(app.GridLayout);
            app.WavelengthLowerLimitSpinner.Step = 20;
            app.WavelengthLowerLimitSpinner.Limits = [200 1500];
            app.WavelengthLowerLimitSpinner.ValueDisplayFormat = '%11.4g nm';
            app.WavelengthLowerLimitSpinner.Layout.Row = 2;
            app.WavelengthLowerLimitSpinner.Layout.Column = 2;
            app.WavelengthLowerLimitSpinner.Value = 450;

            % Create EmissionwavelengthCheckBox
            app.EmissionwavelengthCheckBox = uicheckbox(app.GridLayout);
            app.EmissionwavelengthCheckBox.Text = 'Emission wavelength';
            app.EmissionwavelengthCheckBox.Layout.Row = 2;
            app.EmissionwavelengthCheckBox.Layout.Column = 1;

            % Create SetparameterrangesLabel
            app.SetparameterrangesLabel = uilabel(app.GridLayout);
            app.SetparameterrangesLabel.FontSize = 13;
            app.SetparameterrangesLabel.FontWeight = 'bold';
            app.SetparameterrangesLabel.Layout.Row = 1;
            app.SetparameterrangesLabel.Layout.Column = [1 3];
            app.SetparameterrangesLabel.Text = 'Set parameter ranges';

            % Create NumberofgeneratedPSFsSpinnerLabel
            app.NumberofgeneratedPSFsSpinnerLabel = uilabel(app.GridLayout);
            app.NumberofgeneratedPSFsSpinnerLabel.Layout.Row = 8;
            app.NumberofgeneratedPSFsSpinnerLabel.Layout.Column = 1;
            app.NumberofgeneratedPSFsSpinnerLabel.Text = 'Number of generated PSFs';

            % Create NumberofgeneratedPSFsSpinner
            app.NumberofgeneratedPSFsSpinner = uispinner(app.GridLayout);
            app.NumberofgeneratedPSFsSpinner.Step = 100;
            app.NumberofgeneratedPSFsSpinner.Limits = [1 10000];
            app.NumberofgeneratedPSFsSpinner.RoundFractionalValues = 'on';
            app.NumberofgeneratedPSFsSpinner.ValueDisplayFormat = '%d';
            app.NumberofgeneratedPSFsSpinner.Layout.Row = 8;
            app.NumberofgeneratedPSFsSpinner.Layout.Column = 2;
            app.NumberofgeneratedPSFsSpinner.Value = 2;

            % Create SelectfolderButton
            app.SelectfolderButton = uibutton(app.GridLayout, 'push');
            app.SelectfolderButton.ButtonPushedFcn = createCallbackFcn(app, @SelectfolderButtonPushed, true);
            app.SelectfolderButton.Layout.Row = 9;
            app.SelectfolderButton.Layout.Column = 3;
            app.SelectfolderButton.Text = 'Select folder';

            % Create OutputfolderLabel
            app.OutputfolderLabel = uilabel(app.GridLayout);
            app.OutputfolderLabel.Layout.Row = 9;
            app.OutputfolderLabel.Layout.Column = [1 2];
            app.OutputfolderLabel.Text = 'Output folder:';

            % Create OutputFolder
            app.OutputFolder = uilabel(app.GridLayout);
            app.OutputFolder.Layout.Row = 10;
            app.OutputFolder.Layout.Column = [1 3];
            app.OutputFolder.Text = '';

            % Create SaveimagesCheckBox
            app.SaveimagesCheckBox = uicheckbox(app.GridLayout);
            app.SaveimagesCheckBox.Text = 'Save images';
            app.SaveimagesCheckBox.Layout.Row = 11;
            app.SaveimagesCheckBox.Layout.Column = 1;

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