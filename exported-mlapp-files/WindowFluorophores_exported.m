classdef WindowFluorophores_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        FluorophoresUIFigure  matlab.ui.Figure
        DeleteselectedfluorophoreButton  matlab.ui.control.Button
        AddfluorophoreButton  matlab.ui.control.Button
        FluorophoresLabel     matlab.ui.control.Label
    end

    properties (Access = public)
        TabGroup                     matlab.ui.container.TabGroup
        Tabs                         % cell array of matlab.ui.container.Tab
        xpositionSpinner             % cell array of matlab.ui.control.Spinner
        xpositionSpinnerLabel        % cell array of matlab.ui.control.Label
        ypositionSpinner             % cell array of matlab.ui.control.Spinner
        ypositionSpinnerLabel        % cell array of matlab.ui.control.Label
        zpositionSpinner             % cell array of matlab.ui.control.Spinner
        zpositionSpinnerLabel        % cell array of matlab.ui.control.Label
        InclinationAngleSlider       % cell array of matlab.ui.control.Slider
        InclinationAngleSliderLabel  % cell array of matlab.ui.control.Label
        InclinationAngleEditField    % cell array of matlab.ui.control.EditField
        AzimuthalAngleSlider         % cell array of matlab.ui.control.Slider
        AzimuthalAngleSliderLabel    % cell array of matlab.ui.control.Label
        AzimuthalAngleEditField      % cell array of matlab.ui.control.EditField
        rotationalconstraintSpinner  % cell array of matlab.ui.control.Slider
        rotationalconstraintSpinnerLabel % cell array of matlab.ui.control.Label
        UIAxesPSF                    % cell array of matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp    % Main app object
        psfs          % PSF objects for each fluorophore
        psfImages     % PSF image for each flurophore
    end

    methods (Access = public)

        function n = getNumberFluorophores(app)
            n = numel(app.Tabs);
        end

        function initializePlot(app, psfImage, k)
            % PSF image
            app.psfImages(:,:,k) = psfImage;
            imagesc(app.UIAxesPSF{k}, psfImage);
            axis(app.UIAxesPSF{k}, 'equal');
            axis(app.UIAxesPSF{k}, 'tight');
            cb = colorbar(app.UIAxesPSF{k});
            app.setContrast(k, app.CallingApp.SetcontrastButtonGroup.SelectedObject.Text)
            cb.Label.String = 'Intensity';
            set(app.UIAxesPSF{k},'visible','off')
            colormap(app.UIAxesPSF{k}, app.CallingApp.ColormapDropDown.Value)
        end

        function updatePlot(app, psfImage, k)
            % PSF image
            imagesc(app.UIAxesPSF{k},psfImage);
            app.setContrast(k, app.CallingApp.SetcontrastButtonGroup.SelectedObject.Text)
            colormap(app.UIAxesPSF{k}, app.CallingApp.ColormapDropDown.Value)
        end

        function setContrast(app, k, option)
            switch option
                case 'Intensity limits'
                    psfImage = app.psfImages(:,:,k);
                    upperLimit = max([max(psfImage(:)), 1]);
                    caxis(app.UIAxesPSF{k}, [0 upperLimit])
                case 'Optimize contrast'
                    psfImage = app.psfImages(:,:,k);
                    lowerLimit = max([min(psfImage(:)),0]);
                    upperLimit = max([max(psfImage(:)), lowerLimit + 1]);
                    caxis(app.UIAxesPSF{k}, [lowerLimit upperLimit])
            end
        end

        function setContrastForAll(app, option)
            for k = 1:app.getNumberFluorophores
                switch option
                    case 'Intensity limits'
                        psfImage = app.psfImages(:,:,k);
                        upperLimit = max([max(psfImage(:)), 1]);
                        caxis(app.UIAxesPSF{k}, [0 upperLimit])
                    case 'Optimize contrast'
                        psfImage = app.psfImages(:,:,k);
                        lowerLimit = max([min(psfImage(:)),0]);
                        upperLimit = max([max(psfImage(:)), lowerLimit + 1]);
                        caxis(app.UIAxesPSF{k}, [lowerLimit upperLimit])
                end
            end
        end

        function setColormap(app, map)
            for k = 1:app.getNumberFluorophores
                colormap(app.UIAxesPSF{k}, map)
            end
        end

        function simulateAndDisplayPSFs(app,k)
            app.CallingApp.CalculatingLamp.Color = "y";
            pause(1e-12)
            
            % Get other parameters
            par = app.CallingApp.getParameters();
            
            par.position = Length([app.xpositionSpinner{k}.Value ...
                app.ypositionSpinner{k}.Value ...
                app.zpositionSpinner{k}.Value],'nm'); % position, for xy position 0 corresponds to the center of the center pixel

            theta = app.InclinationAngleSlider{k}.Value; % dipole inclination angle (in degree)
            phi = app.AzimuthalAngleSlider{k}.Value; % dipole azimuthal angle
            par.dipole = Dipole(theta*pi/180, phi*pi/ 180); % conversion to rad
            par.rotationalConstraint = app.rotationalconstraintSpinner{k}.Value;
            

            % Simulate PSFs of fluorophores
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

            app.psfs{k} = psf;

            if size(app.psfImages,1) == size(psf.image,1)
                app.psfImages(:,:,k) = psf.image;
            else
                app.psfImages = [];
                app.psfImages(:,:,k) = psf.image;
            end
            app.updatePlot(app.psfImages(:,:,k),k)
            
            if app.CallingApp.ShowPSFCheckBox.Value
                app.updateCombinedPSF()
            end
            if app.CallingApp.PolarizedemissionchannelsCheckBox.Value
                app.updatePolarizedEmissionChannels()
            end
            app.CallingApp.CalculatingLamp.Color = "g";
        end

        function updateCombinedPSF(app)
            sumImage = sum(app.psfImages,3);
            if app.CallingApp.ShowPsf2DCheckBox.Value
                app.CallingApp.PlotPSF.updatePlot(sumImage);
            end
        end

        function updatePolarizedEmissionChannels(app)
            Ixs = cellfun(@(x) x.Ix, app.psfs, 'UniformOutput', false);
            psf.Ix = sum(cat(3,Ixs{:}),3);
            Iys = cellfun(@(x) x.Iy, app.psfs, 'UniformOutput', false);
            psf.Iy = sum(cat(3,Iys{:}),3);
            app.CallingApp.PlotPolarizedEmission.updatePlot(psf);
        end

        function updateDipoleRotation(app,str,k)
            if strcmp(str, 'fixed')
                app.AzimuthalAngleSliderLabel{k}.Visible = 'on';
                app.AzimuthalAngleSlider{k}.Visible = 'on';
                app.AzimuthalAngleEditField{k}.Visible = 'on';
                app.InclinationAngleSliderLabel{k}.Visible = 'on';
                app.InclinationAngleSlider{k}.Visible = 'on';
                app.InclinationAngleEditField{k}.Visible = 'on';
                app.rotationalconstraintSpinner{k}.Visible = 'off';
                app.rotationalconstraintSpinnerLabel{k}.Visible = 'off';
            elseif strcmp(str,'freely rotating')
                app.AzimuthalAngleSliderLabel{k}.Visible = 'off';
                app.AzimuthalAngleSlider{k}.Visible = 'off';
                app.AzimuthalAngleEditField{k}.Visible = 'off';
                app.InclinationAngleSliderLabel{k}.Visible = 'off';
                app.InclinationAngleSlider{k}.Visible = 'off';
                app.InclinationAngleEditField{k}.Visible = 'off';
                app.rotationalconstraintSpinner{k}.Visible = 'off';
                app.rotationalconstraintSpinnerLabel{k}.Visible = 'off';
            elseif strcmp(str, 'partially rotating')
                app.AzimuthalAngleSliderLabel{k}.Visible = 'on';
                app.AzimuthalAngleSlider{k}.Visible = 'on';
                app.AzimuthalAngleEditField{k}.Visible = 'on';
                app.InclinationAngleSliderLabel{k}.Visible = 'on';
                app.InclinationAngleSlider{k}.Visible = 'on';
                app.InclinationAngleEditField{k}.Visible = 'on';
                app.rotationalconstraintSpinner{k}.Visible = 'on';
                app.rotationalconstraintSpinnerLabel{k}.Visible = 'on';
            else
                error('Dipole rotation must be either fixed, freely rotating or partially rotating!')
            end
        end
    end


    methods (Access=private)
        function createNewFluorophoreTab(app,k)
            app.Tabs{k} = uitab(app.TabGroup);
            app.Tabs{k}.Title = num2str(k);

            %% PSF image
            % Create UIAxesPSF
            app.UIAxesPSF{k} = uiaxes(app.Tabs{k});
            app.UIAxesPSF{k}.PlotBoxAspectRatio = [1 1 1];
            app.UIAxesPSF{k}.XTick = [];
            app.UIAxesPSF{k}.XTickLabelRotation = 0;
            app.UIAxesPSF{k}.YTick = [];
            app.UIAxesPSF{k}.YTickLabelRotation = 0;
            app.UIAxesPSF{k}.ZTickLabelRotation = 0;
            app.UIAxesPSF{k}.FontSize = 15;
            app.UIAxesPSF{k}.Position = [345 25 280 231];

            %% Position

            % Create xpositionSpinnerLabel
            app.xpositionSpinnerLabel{k} = uilabel(app.Tabs{k});
            app.xpositionSpinnerLabel{k}.Position = [17 233 57 22];
            
            app.xpositionSpinnerLabel{k}.Text = 'x-position';
            app.xpositionSpinnerLabel{k}.BusyAction = 'cancel';

            % Create xpositionSpinner
            app.xpositionSpinner{k} = uispinner(app.Tabs{k});
            app.xpositionSpinner{k}.Limits = [-inf inf];
            app.xpositionSpinner{k}.Step = 10;
            app.xpositionSpinner{k}.ValueDisplayFormat = '%11.4g nm';
            app.xpositionSpinner{k}.ValueChangedFcn = {@app.xpositionSpinnerValueChanged, k};
            app.xpositionSpinner{k}.ValueChangingFcn = {@app.xpositionSpinnerValueChanging, k};
            app.xpositionSpinner{k}.Position = [130 233 102 22];
            app.xpositionSpinner{k}.BusyAction = 'cancel';

            % Create ypositionSpinnerLabel
            app.ypositionSpinnerLabel{k} = uilabel(app.Tabs{k});
            app.ypositionSpinnerLabel{k}.Position = [17 203 57 22];
            app.ypositionSpinnerLabel{k}.Text = 'y-position';
            app.ypositionSpinnerLabel{k}.BusyAction = 'cancel';

            % Create ypositionSpinner
            app.ypositionSpinner{k} = uispinner(app.Tabs{k});
            app.ypositionSpinner{k}.Limits = [-inf inf];
            app.ypositionSpinner{k}.Step = 10;
            app.ypositionSpinner{k}.ValueDisplayFormat = '%11.4g nm';
            app.ypositionSpinner{k}.ValueChangedFcn = {@app.ypositionSpinnerValueChanged, k};
            app.ypositionSpinner{k}.ValueChangingFcn = {@app.ypositionSpinnerValueChanging, k};
            app.ypositionSpinner{k}.Position = [130 203 102 22];
            app.ypositionSpinner{k}.BusyAction = 'cancel';

            % Create zpositionSpinnerLabel
            app.zpositionSpinnerLabel{k} = uilabel(app.Tabs{k});
            app.zpositionSpinnerLabel{k}.Position = [17 173 57 22];
            app.zpositionSpinnerLabel{k}.Text = 'z-position';
            app.zpositionSpinnerLabel{k}.BusyAction = 'cancel';

            % Create zpositionSpinner
            app.zpositionSpinner{k} = uispinner(app.Tabs{k});
            app.zpositionSpinner{k}.Limits = [0 inf];
            app.zpositionSpinner{k}.Step = 10;
            app.zpositionSpinner{k}.ValueDisplayFormat = '%11.4g nm';
            app.zpositionSpinner{k}.ValueChangedFcn = {@app.zpositionSpinnerValueChanged, k};
            app.zpositionSpinner{k}.ValueChangingFcn = {@app.zpositionSpinnerValueChanging, k};
            app.zpositionSpinner{k}.Position = [130 173 102 22];
            app.zpositionSpinner{k}.BusyAction = 'cancel';
           
            %% Dipole orientation

            % Create AzimuthalAngleSliderLabel
            app.AzimuthalAngleSliderLabel{k} = uilabel(app.Tabs{k});
            app.AzimuthalAngleSliderLabel{k}.Position = [17 73 92 22];
            app.AzimuthalAngleSliderLabel{k}.Text = 'Azimuthal Angle';
            app.AzimuthalAngleSliderLabel{k}.BusyAction = 'cancel';

            % Create AzimuthalAngleSlider
            app.AzimuthalAngleSlider{k} = uislider(app.Tabs{k});
            app.AzimuthalAngleSlider{k}.Limits = [0 360];
            app.AzimuthalAngleSlider{k}.ValueChangedFcn = {@app.AzimuthalAngleSliderValueChanged, k};
            app.AzimuthalAngleSlider{k}.ValueChangingFcn = {@app.AzimuthalAngleSliderValueChanging, k};
            app.AzimuthalAngleSlider{k}.Position = [130 82 140 3];
            app.AzimuthalAngleSlider{k}.BusyAction = 'cancel';

            % Create AzimuthalAngleEditField
            app.AzimuthalAngleEditField{k} = uieditfield(app.Tabs{k}, 'numeric');
            app.AzimuthalAngleEditField{k}.Limits = [0 360];
            app.AzimuthalAngleEditField{k}.RoundFractionalValues = 'on';
            app.AzimuthalAngleEditField{k}.ValueDisplayFormat = '%d°';
            app.AzimuthalAngleEditField{k}.ValueChangedFcn = {@app.AzimuthalAngleEditFieldValueChanged, k};
            app.AzimuthalAngleEditField{k}.Position = [290 73 43 22];

            % Create InclinationAngleSliderLabel
            app.InclinationAngleSliderLabel{k} = uilabel(app.Tabs{k});
            app.InclinationAngleSliderLabel{k}.Position = [17 125 94 22];
            app.InclinationAngleSliderLabel{k}.Text = 'Inclination Angle';
            app.InclinationAngleSliderLabel{k}.BusyAction = 'cancel';

            % Create InclinationAngleSlider
            app.InclinationAngleSlider{k} = uislider(app.Tabs{k});
            app.InclinationAngleSlider{k}.Limits = [0 90];
            app.InclinationAngleSlider{k}.ValueChangedFcn = {@app.InclinationAngleSliderValueChanged, k};
            app.InclinationAngleSlider{k}.ValueChangingFcn = {@app.InclinationAngleSliderValueChanging, k};
            app.InclinationAngleSlider{k}.Position = [130 134 140 3];
            app.InclinationAngleSlider{k}.BusyAction = 'cancel';
            
            % Create InclinationAngleEditField
            app.InclinationAngleEditField{k} = uieditfield(app.Tabs{k}, 'numeric');
            app.InclinationAngleEditField{k}.Limits = [0 90];
            app.InclinationAngleEditField{k}.RoundFractionalValues = 'on';
            app.InclinationAngleEditField{k}.ValueDisplayFormat = '%d°';
            app.InclinationAngleEditField{k}.ValueChangedFcn = {@app.InclinationAngleEditFieldValueChanged, k};
            app.InclinationAngleEditField{k}.Position = [290 125 43 22];
            

            % Create rotationalconstraintSpinnerLabel
            app.rotationalconstraintSpinnerLabel{k} = uilabel(app.Tabs{k});
            app.rotationalconstraintSpinnerLabel{k}.Position = [17 10 130 22];
            app.rotationalconstraintSpinnerLabel{k}.Text = 'Rotational constraint';
            app.rotationalconstraintSpinnerLabel{k}.BusyAction = 'cancel';

            % Create rotationalconstraintSpinner
            app.rotationalconstraintSpinner{k} = uispinner(app.Tabs{k});
            app.rotationalconstraintSpinner{k}.Limits = [0 1];
            app.rotationalconstraintSpinner{k}.Step = 0.05;
            app.rotationalconstraintSpinner{k}.Value = 0.1;
            app.rotationalconstraintSpinner{k}.ValueDisplayFormat = '%11.4g';
            app.rotationalconstraintSpinner{k}.ValueChangedFcn = {@app.rotationalconstraintSpinnerValueChanged, k};
            app.rotationalconstraintSpinner{k}.ValueChangingFcn = {@app.rotationalconstraintSpinnerValueChanging, k};
            app.rotationalconstraintSpinner{k}.Position = [150 10 60 22];
            app.rotationalconstraintSpinner{k}.BusyAction = 'cancel';
            
            % Hide spinners if freely rotating dipole is selected
            app.updateDipoleRotation(app.CallingApp.DipolerotationButtonGroup.SelectedObject.Text,k)
            
            % Hide delete button if there is only one fluorophore (at beginning)
            if app.getNumberFluorophores==1
                app.DeleteselectedfluorophoreButton.Visible = "off";
            end
        end

        function DeletefluorophoreButtonPushed(app, k)
            
            % Delete PSF image
            app.psfImages(:,:,k) = [];

            % Recalculate PSF sum of remaining fluorophores
            app.updateCombinedPSF()
            
            % Delete tab
            delete(app.Tabs{k});
            app.Tabs{k} = []; % make sure nowhere else will refer to a deleted object
            app.Tabs(k) = []; % app.Tabs(remainingTabs); % remove from cell array
            
            delete(app.xpositionSpinner{k});
            app.xpositionSpinner{k} = [];  
            app.xpositionSpinner(k) = []; 

            delete(app.xpositionSpinnerLabel{k});
            app.xpositionSpinnerLabel{k} = [];  
            app.xpositionSpinnerLabel(k) = []; 

            delete(app.ypositionSpinner{k});
            app.ypositionSpinner{k} = [];  
            app.ypositionSpinner(k) = []; 

            delete(app.ypositionSpinnerLabel{k});
            app.ypositionSpinnerLabel{k} = [];  
            app.ypositionSpinnerLabel(k) = [];

            delete(app.zpositionSpinner{k});
            app.zpositionSpinner{k} = [];  
            app.zpositionSpinner(k) = [];

            delete(app.zpositionSpinnerLabel{k});
            app.zpositionSpinnerLabel{k} = [];  
            app.zpositionSpinnerLabel(k) = [];
            
            delete(app.InclinationAngleSlider{k});
            app.InclinationAngleSlider{k} = [];  
            app.InclinationAngleSlider(k) = []; 

            delete(app.InclinationAngleSliderLabel{k});
            app.InclinationAngleSliderLabel{k} = [];  
            app.InclinationAngleSliderLabel(k) = []; 

            delete(app.AzimuthalAngleSlider{k});
            app.AzimuthalAngleSlider{k} = [];  
            app.AzimuthalAngleSlider(k) = []; 

            delete(app.AzimuthalAngleSliderLabel{k});
            app.AzimuthalAngleSliderLabel{k} = [];  
            app.AzimuthalAngleSliderLabel(k) = []; 
            
            delete(app.InclinationAngleEditField{k});
            app.InclinationAngleEditField{k} = [];  
            app.InclinationAngleEditField(k) = [];

            delete(app.AzimuthalAngleEditField{k});
            app.AzimuthalAngleEditField{k} = [];  
            app.AzimuthalAngleEditField(k) = [];
            
            delete(app.rotationalconstraintSpinner{k});
            app.rotationalconstraintSpinner{k} = [];  
            app.rotationalconstraintSpinner(k) = []; 

            delete(app.rotationalconstraintSpinnerLabel{k});
            app.rotationalconstraintSpinnerLabel{k} = [];  
            app.rotationalconstraintSpinnerLabel(k) = [];
            
            delete(app.UIAxesPSF{k});
            app.UIAxesPSF{k} = [];  
            app.UIAxesPSF(k) = [];

            % update function handles
            for index = k:app.getNumberFluorophores
                app.Tabs{index}.Title = num2str(index);
                app.xpositionSpinner{index}.ValueChangedFcn = {@app.xpositionSpinnerValueChanged, index};
                app.xpositionSpinner{index}.ValueChangingFcn = {@app.xpositionSpinnerValueChanging, index};
                app.ypositionSpinner{index}.ValueChangedFcn = {@app.ypositionSpinnerValueChanged, index};
                app.ypositionSpinner{index}.ValueChangingFcn = {@app.ypositionSpinnerValueChanging, index};
                app.zpositionSpinner{index}.ValueChangedFcn = {@app.zpositionSpinnerValueChanged, index};
                app.zpositionSpinner{index}.ValueChangingFcn = {@app.zpositionSpinnerValueChanging, index};
                app.AzimuthalAngleSlider{index}.ValueChangedFcn = {@app.AzimuthalAngleSliderValueChanged, index};
                app.AzimuthalAngleSlider{index}.ValueChangingFcn = {@app.AzimuthalAngleSliderValueChanging, index};
                app.AzimuthalAngleEditField{index}.ValueChangedFcn = {@app.AzimuthalAngleEditFieldValueChanged, index};
                app.InclinationAngleSlider{index}.ValueChangedFcn = {@app.InclinationAngleSliderValueChanged, index};
                app.InclinationAngleSlider{index}.ValueChangingFcn = {@app.InclinationAngleSliderValueChanging, index};
                app.InclinationAngleEditField{index}.ValueChangedFcn = {@app.InclinationAngleEditFieldValueChanged, index};
                app.rotationalconstraintSpinner{index}.ValueChangedFcn = {@app.rotationalconstraintSpinnerValueChanged, index};
                app.rotationalconstraintSpinner{index}.ValueChangingFcn = {@app.rotationalconstraintSpinnerValueChanging, index};
            end

            % Make delete button invisible if there is only one fluorophore remaining 
            if app.getNumberFluorophores==1
                app.DeleteselectedfluorophoreButton.Visible = "off";
            end
        end

        function xpositionSpinnerValueChanged(app, src, event, k)
            app.xpositionSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function xpositionSpinnerValueChanging(app, src, event, k)
            app.xpositionSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function ypositionSpinnerValueChanged(app, src, event, k)
            app.ypositionSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function ypositionSpinnerValueChanging(app, src, event, k)
            app.ypositionSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function zpositionSpinnerValueChanged(app, src, event, k)
            app.zpositionSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function zpositionSpinnerValueChanging(app, src, event, k)
            app.zpositionSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function InclinationAngleSliderValueChanged(app, src, event, k)
            app.InclinationAngleSlider{k}.Value = event.Value;
            app.InclinationAngleEditField{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function InclinationAngleSliderValueChanging(app, src, event, k)
            app.InclinationAngleSlider{k}.Value = event.Value;
            app.InclinationAngleEditField{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function InclinationAngleEditFieldValueChanged(app, src, event, k)
            app.InclinationAngleEditField{k}.Value = event.Value;
            app.InclinationAngleSlider{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k);
        end

        function AzimuthalAngleSliderValueChanged(app, src, event, k)
            app.AzimuthalAngleSlider{k}.Value = event.Value;
            app.AzimuthalAngleEditField{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function AzimuthalAngleSliderValueChanging(app, src, event, k)
            app.AzimuthalAngleSlider{k}.Value = event.Value;
            app.AzimuthalAngleEditField{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function AzimuthalAngleEditFieldValueChanged(app, src, event, k)
            app.AzimuthalAngleEditField{k}.Value = event.Value;
            app.AzimuthalAngleSlider{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k);
        end
        
        function rotationalconstraintSpinnerValueChanged(app, src, event, k)
            app.rotationalconstraintSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end

        function rotationalconstraintSpinnerValueChanging(app, src, event, k)
            app.rotationalconstraintSpinner{k}.Value = event.Value;
            simulateAndDisplayPSFs(app,k)
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            app.CallingApp.CalculatingLamp.Color = "y";
            pause(1e-12)
            
            % Store main app object
            app.CallingApp = mainapp;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.FluorophoresUIFigure);
            app.TabGroup.Position = [17 17 656 295];

            % Create Tabs for first 2 fluorophores
            app.Tabs = cell(1,1);
            for k = 1:1
                app.createNewFluorophoreTab(k)
                app.initializePlot(zeros(app.CallingApp.nPixels), k)
                app.simulateAndDisplayPSFs(k)
            end

            app.CallingApp.CalculatingLamp.Color = "g";
        end

        % Close request function: FluorophoresUIFigure
        function FluorophoresUIFigureCloseRequest(app, event)
            app.CallingApp.SwitchMultipleFluorophores = 0;
            app.CallingApp.AzimuthalAngleSliderLabel.Visible = "on";
            app.CallingApp.AzimuthalAngleSlider.Visible = "on";
            app.CallingApp.phi.Visible = "on";
            app.CallingApp.InclinationAngleSliderLabel.Visible = "on";
            app.CallingApp.InclinationAngleSlider.Visible = "on";
            app.CallingApp.theta.Visible = "on";
            app.CallingApp.xpositionSpinnerLabel.Visible = "on";
            app.CallingApp.xpositionSpinner.Visible = "on";
            app.CallingApp.ypositionSpinnerLabel.Visible = "on";
            app.CallingApp.ypositionSpinner.Visible = "on";
            app.CallingApp.zpositionSpinnerLabel.Visible = "on";
            app.CallingApp.zpositionSpinner.Visible = "on";
            app.CallingApp.PositionLabel.Visible = 'on';
            app.CallingApp.DipolerotationButtonGroup.Visible = 'on';
            app.CallingApp.DipoleorientationLabel.Visible = 'on';
            app.CallingApp.simulateAndDisplayPSF();
            set(app.CallingApp.ConfiguremultiplefluorophoresButton, 'Enable', 'on')

            app.CallingApp.CalculateCheckBox.Visible = 'on'; % checkbox for CRB
            app.CallingApp.CRBOutputField.Visible = 'on';
            app.CallingApp.CramrRaoBoundLabel.Visible = 'on';
            % Switch on option for 3D plot again (currently not compatible
            % with multiple fluorophores
            if app.CallingApp.ShowPSFCheckBox.Value == 1
                set(app.CallingApp.ShowPsf3DCheckBox, 'Enable', 'on')
            end
            delete(app)
        end

        % Button pushed function: AddfluorophoreButton
        function AddfluorophoreButtonPushed(app, event)
            % Get number of new fluorophore
            k = app.getNumberFluorophores + 1;

            app.createNewFluorophoreTab(k)
            app.initializePlot(zeros(app.CallingApp.nPixels), k)
            app.simulateAndDisplayPSFs(k)
            
            app.TabGroup.SelectedTab = app.Tabs{k};

            % enable deletion of first fluorophore 
            if app.getNumberFluorophores>1
                app.DeleteselectedfluorophoreButton.Visible = "on";
            end
        end

        % Button pushed function: DeleteselectedfluorophoreButton
        function DeleteselectedfluorophoreButtonPushed(app, event)
            % Get number of currently selected fluorophore
            k = int32(str2double(app.TabGroup.SelectedTab.Title));
            app.DeletefluorophoreButtonPushed(k)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create FluorophoresUIFigure and hide until all components are created
            app.FluorophoresUIFigure = uifigure('Visible', 'off');
            app.FluorophoresUIFigure.Position = [550 80 691 372];
            app.FluorophoresUIFigure.Name = 'Fluorophores';
            app.FluorophoresUIFigure.CloseRequestFcn = createCallbackFcn(app, @FluorophoresUIFigureCloseRequest, true);
            app.FluorophoresUIFigure.Scrollable = 'on';

            % Create FluorophoresLabel
            app.FluorophoresLabel = uilabel(app.FluorophoresUIFigure);
            app.FluorophoresLabel.FontSize = 13;
            app.FluorophoresLabel.FontWeight = 'bold';
            app.FluorophoresLabel.Position = [19 333 89 22];
            app.FluorophoresLabel.Text = 'Fluorophores';

            % Create AddfluorophoreButton
            app.AddfluorophoreButton = uibutton(app.FluorophoresUIFigure, 'push');
            app.AddfluorophoreButton.ButtonPushedFcn = createCallbackFcn(app, @AddfluorophoreButtonPushed, true);
            app.AddfluorophoreButton.Position = [573 333 101 22];
            app.AddfluorophoreButton.Text = 'Add fluorophore';

            % Create DeleteselectedfluorophoreButton
            app.DeleteselectedfluorophoreButton = uibutton(app.FluorophoresUIFigure, 'push');
            app.DeleteselectedfluorophoreButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteselectedfluorophoreButtonPushed, true);
            app.DeleteselectedfluorophoreButton.Position = [400 333 162 23];
            app.DeleteselectedfluorophoreButton.Text = 'Delete selected fluorophore';

            % Show the figure after all components are created
            app.FluorophoresUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = WindowFluorophores_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.FluorophoresUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.FluorophoresUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.FluorophoresUIFigure)
        end
    end
end