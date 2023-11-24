classdef MainSimulationPSF_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PSFsimulationUIFigure           matlab.ui.Figure
        CalculatingLamp                 matlab.ui.control.Lamp
        TabGroup                        matlab.ui.container.TabGroup
        FluorophoreTab                  matlab.ui.container.Tab
        OrientationLabel                matlab.ui.control.Label
        PositionLabel                   matlab.ui.control.Label
        xpositionSpinner                matlab.ui.control.Spinner
        xpositionSpinnerLabel           matlab.ui.control.Label
        ypositionSpinner                matlab.ui.control.Spinner
        ypositionSpinnerLabel           matlab.ui.control.Label
        EmissionWavelengthSpinner       matlab.ui.control.Spinner
        EmissionwavelengthSpinnerLabel  matlab.ui.control.Label
        DipoleorientationLabel          matlab.ui.control.Label
        DipolerotationButtonGroup       matlab.ui.container.ButtonGroup
        fixedButton                     matlab.ui.control.RadioButton
        freelyrotatingButton            matlab.ui.control.RadioButton
        phi                             matlab.ui.control.NumericEditField
        theta                           matlab.ui.control.NumericEditField
        zpositionSpinner                matlab.ui.control.Spinner
        zpositionSpinnerLabel           matlab.ui.control.Label
        ConfiguremultiplefluorophoresButton  matlab.ui.control.Button
        ReducedexcitationSwitch         matlab.ui.control.Switch
        ReducedexcitationSwitchLabel    matlab.ui.control.Label
        PhotonnumberSpinner             matlab.ui.control.Spinner
        PhotonnumberSpinnerLabel        matlab.ui.control.Label
        PhotonshotnoiseSwitch           matlab.ui.control.Switch
        PhotonshotnoiseSwitchLabel      matlab.ui.control.Label
        InclinationAngleSlider          matlab.ui.control.Slider
        InclinationAngleSliderLabel     matlab.ui.control.Label
        AzimuthalAngleSlider            matlab.ui.control.Slider
        AzimuthalAngleSliderLabel       matlab.ui.control.Label
        FluorophoreLabel                matlab.ui.control.Label
        MicroscopeRITab                 matlab.ui.container.Tab
        ObjectiveFocalLengthSpinner     matlab.ui.control.Spinner
        ObjectivefocallengthLabel       matlab.ui.control.Label
        TubeLensFocalLengthSpinner      matlab.ui.control.Spinner
        TubelensfocallengthLabel        matlab.ui.control.Label
        PixelsizeObjectSpaceSpinner     matlab.ui.control.Spinner
        PixelsizeobjectspaceLabel       matlab.ui.control.Label
        PixelsizePhysicalSpinner        matlab.ui.control.Spinner
        PixelsizephysicalSpinnerLabel   matlab.ui.control.Label
        BackgroundnoisestdSpinner       matlab.ui.control.Spinner
        BackgroundnoisestdSpinnerLabel  matlab.ui.control.Label
        CameraLabel                     matlab.ui.control.Label
        IntermediateLayerThicknessSpinner  matlab.ui.control.Spinner
        IntermediateLayerThicknessLabel  matlab.ui.control.Label
        AddintermediatelayerCheckBox    matlab.ui.control.CheckBox
        MagnificationSpinner            matlab.ui.control.Spinner
        MagnificationLabel              matlab.ui.control.Label
        RefractiveIndexImmersionMediumSpinner  matlab.ui.control.Spinner
        RefractiveIndexImmersionMediumLabel  matlab.ui.control.Label
        RefractiveIndexIntermediateLayerSpinner  matlab.ui.control.Spinner
        RefractiveIndexIntermediateLayerLabel  matlab.ui.control.Label
        RefractiveIndexSpecimenSpinner  matlab.ui.control.Spinner
        RefractiveIndexSpecimenLabel    matlab.ui.control.Label
        RefractiveindexLabel            matlab.ui.control.Label
        defocus                         matlab.ui.control.NumericEditField
        FocusButton                     matlab.ui.control.Button
        ObjectiveNaSpinner              matlab.ui.control.Spinner
        ObjectiveNALabel                matlab.ui.control.Label
        DefocusSlider                   matlab.ui.control.Slider
        DefocusSliderLabel              matlab.ui.control.Label
        ObjectiveNaLabel                matlab.ui.control.Label
        ObjectiveandtubelensLabel       matlab.ui.control.Label
        AberrationsTab                  matlab.ui.container.Tab
        WeightHorizontalComaSpinner     matlab.ui.control.Spinner
        WeightHorizontalComaSpinnerLabel  matlab.ui.control.Label
        WeightVerticalComaSpinner       matlab.ui.control.Spinner
        WeightVerticalComaSpinnerLabel  matlab.ui.control.Label
        WeightVerticalAstigmatismSpinner  matlab.ui.control.Spinner
        WeightVerticalAstigmatismSpinnerLabel  matlab.ui.control.Label
        WeightObliqueAstigmatismSpinner  matlab.ui.control.Spinner
        WeightObliqueAstigmatismSpinnerLabel  matlab.ui.control.Label
        WeightPrimarySphericalSpinner   matlab.ui.control.Spinner
        WeightPrimarySphericalSpinnerLabel  matlab.ui.control.Label
        WeightDefocusSpinner            matlab.ui.control.Spinner
        WeightDefocusSpinnerLabel       matlab.ui.control.Label
        WeightHorizontalTiltSpinner     matlab.ui.control.Spinner
        WeightHorizontalTiltSpinnerLabel  matlab.ui.control.Label
        VerticalTiltCheckBox            matlab.ui.control.CheckBox
        WeightVerticalTiltSpinner       matlab.ui.control.Spinner
        WeightVerticalTiltSpinnerLabel  matlab.ui.control.Label
        HorizontalTiltCheckBox          matlab.ui.control.CheckBox
        DefocusCheckBox                 matlab.ui.control.CheckBox
        PrimarySphericalCheckBox        matlab.ui.control.CheckBox
        ObliqueAstigmatismCheckBox      matlab.ui.control.CheckBox
        VerticalAstigmatismCheckBox     matlab.ui.control.CheckBox
        VerticalComaCheckBox            matlab.ui.control.CheckBox
        HorizontalComaCheckBox          matlab.ui.control.CheckBox
        FitZernikeButton                matlab.ui.control.Button
        FitZernikeLabel                 matlab.ui.control.Label
        FitZernikeaberrationsLabel      matlab.ui.control.Label
        ZernikeAberrationsShowplotCheckBox  matlab.ui.control.CheckBox
        ZernikeaberrationsLabel         matlab.ui.control.Label
        SpecifyZernikeaberrationsButtonGroup  matlab.ui.container.ButtonGroup
        SaveZernikeIndicesWeightsButton  matlab.ui.control.Button
        IndicesWeightsEditField         matlab.ui.control.EditField
        LoadZernikeIndicesWeightsButton  matlab.ui.control.Button
        ZernikeInputVectorButton        matlab.ui.control.RadioButton
        ZernikeCommonAberrationsButton  matlab.ui.control.RadioButton
        PhasemaskTab                    matlab.ui.container.Tab
        PhaseMaskFilepathLabel          matlab.ui.control.Label
        LoadCustomPhaseMaskButton       matlab.ui.control.Button
        PhaseMaskShowplotCheckBox       matlab.ui.control.CheckBox
        MaxShiftSlider                  matlab.ui.control.Slider
        MaxShiftSliderLabel             matlab.ui.control.Label
        NumberfacetsSpinner             matlab.ui.control.Spinner
        NumberfacetsSpinnerLabel        matlab.ui.control.Label
        PhaseMaskOptionsLabel           matlab.ui.control.Label
        SectorSlider                    matlab.ui.control.Slider
        SectorSliderLabel               matlab.ui.control.Label
        RotatephasemaskSlider           matlab.ui.control.Slider
        RotatephasemaskSliderLabel      matlab.ui.control.Label
        InnerringradiusSlider           matlab.ui.control.Slider
        InnerringradiusSliderLabel      matlab.ui.control.Label
        BFPmanipulationDropDown         matlab.ui.control.DropDown
        BFPmanipulationDropDownLabel    matlab.ui.control.Label
        PhasemaskLabel                  matlab.ui.control.Label
        TransmissionTab                 matlab.ui.container.Tab
        TransmissionMaskFilepathLabel   matlab.ui.control.Label
        LoadCustomTransmissionMaskButton  matlab.ui.control.Button
        TransmissionmaskDropDown        matlab.ui.control.DropDown
        TransmissionmaskDropDownLabel   matlab.ui.control.Label
        TransmissionMaskShowplotCheckBox  matlab.ui.control.CheckBox
        TransmissionLabel               matlab.ui.control.Label
        OptionsTab                      matlab.ui.container.Tab
        CramrRaoBoundLabel              matlab.ui.control.Label
        CRBOutputField                  matlab.ui.control.Label
        CalculateCramrRaoBoundCheckBox  matlab.ui.control.CheckBox
        StepsizeEditField               matlab.ui.control.NumericEditField
        StepsizeEditFieldLabel          matlab.ui.control.Label
        NumberstepsEditField            matlab.ui.control.NumericEditField
        NumberstepsEditFieldLabel       matlab.ui.control.Label
        ShowPsf2DCheckBox               matlab.ui.control.CheckBox
        ShowPsf3DCheckBox               matlab.ui.control.CheckBox
        PixelsperlateralaxisEditField   matlab.ui.control.NumericEditField
        PixelsperlateralaxisEditFieldLabel  matlab.ui.control.Label
        ContrastLabel                   matlab.ui.control.Label
        PlotoptionsLabel                matlab.ui.control.Label
        WindowselectionLabel            matlab.ui.control.Label
        ShowPSFCheckBox                 matlab.ui.control.CheckBox
        PolarizedemissionchannelsCheckBox  matlab.ui.control.CheckBox
        ROIsidelengthEditField          matlab.ui.control.NumericEditField
        ROIsidelengthEditFieldLabel     matlab.ui.control.Label
        ColormapDropDown                matlab.ui.control.DropDown
        ColormapDropDownLabel           matlab.ui.control.Label
        SetcontrastButtonGroup          matlab.ui.container.ButtonGroup
        OptimizecontrastButton          matlab.ui.control.RadioButton
        IntensitylimitsButton           matlab.ui.control.RadioButton
        SetparametersLabel              matlab.ui.control.Label
    end

    properties (Access = public)
        PlotPSF             % Plot PSF app
        PlotPSFThreeDim             % Plot PSF 3D app
        PlotPhaseMask       % Plot Phase mask app
        PlotTransmissionMask
        PlotZernikeAberrations % Plot Zernike aberrations app
        PlotPolarizedEmission % Plot Polarized emission app
        Fluorophores        % Fluorophore window
        SwitchMultipleFluorophores logical = 0
        FitZernikeApp % Call app for Zernike fitting

        phaseMask % store loaded phase mask
        transmissionMask % store loaded transmission mask
        zernikeIndices % store Zernike indices
        zernikeWeights % store Zernike weights

        nPixels {mustBeInteger} = 25
        nDiscretizationBFP = 129
        
        parameters
    end
    
    properties (Hidden)
        originalPath
        isTransmissionMaskLoaded = false
    end

    methods (Access = public)
        function psf = simulateAndDisplayPSF(app)

            app.CalculatingLamp.Color = "y";
            pause(1e-12)
            
            par = getParameters(app);
            app.parameters = par;
            
            % Update phase mask, transmission and Zernike aberration plots if checkbox is activated 
            if app.PhaseMaskShowplotCheckBox.Value
                app.PlotPhaseMask.updatePlot();
            end
            if app.TransmissionMaskShowplotCheckBox.Value
                app.PlotTransmissionMask.updatePlot();
            end
            if app.ZernikeAberrationsShowplotCheckBox.Value
                app.PlotZernikeAberrations.updatePlot();
            end

            % Flurophores
            switch app.SwitchMultipleFluorophores
                case 0 % 'Single'
                    par.dipole = Dipole(app.theta.Value*pi/180, app.phi.Value*pi/ 180); % conversion to rad, theta = inclination angle, phi = azimuthal angle in degree
                    par.position = Length([app.xpositionSpinner.Value app.ypositionSpinner.Value app.zpositionSpinner.Value],'nm'); % position, for xy position 0 corresponds to the center of the center pixel

                    % Simulate PSF
                    switch app.DipolerotationButtonGroup.SelectedObject.Text
                        case 'fixed'
                            psf = PSF(par);
                        case 'freely rotating'
                            psf = IsotropicPSF(par);
                        otherwise
                            error('Invalid input value for dipole rotation!')
                    end

                    psfImage = psf.image;
                    if app.ShowPSFCheckBox.Value
                        if numel(par.defocus)>1
                            nSteps = numel(par.defocus);
                            midSlice = ceil(nSteps/2);
                            app.PlotPSF.updatePlot(psfImage(:,:,midSlice));
                        else
                            app.PlotPSF.updatePlot(psfImage);
                        end
                    end
                    if app.ShowPsf3DCheckBox.Value
                        app.PlotPSFThreeDim.updatePlot(psfImage);
                    end
                    % update polarized channels plot if checkbox is activated 
                    if app.PolarizedemissionchannelsCheckBox.Value
                        app.PlotPolarizedEmission.updatePlot(psf);
                    end

                case 1 % 'Multiple'
                    for k = 1:app.Fluorophores.getNumberFluorophores
                        app.Fluorophores.simulateAndDisplayPSFs(k)
                    end
            end
            
            if strcmp(app.TabGroup.SelectedTab.Title, 'Options') && app.CalculateCramrRaoBoundCheckBox.Value 
                psf.CRB = psf.calculateCRB;
                CRBinNanometer = round(psf.CRB.inNanometer, 3, 'significant');
                upperBound = 200; % display upper bound if CRB is larger than some threshold
                if CRBinNanometer(1) > upperBound
                    CRBxOutput = ['x:  ','> ', num2str(200), ' nm'];
                else
                    CRBxOutput = ['x:  ',num2str(CRBinNanometer(1)), ' nm']; 
                end

                if CRBinNanometer(2) > upperBound
                    CRByOutput = ['y:  ','> ', num2str(200), ' nm'];
                else
                    CRByOutput = ['y:  ',num2str(CRBinNanometer(2)), ' nm']; 
                end

                if CRBinNanometer(3) > upperBound
                    CRBzOutput = ['z:  ','> ', num2str(200), ' nm'];
                else
                    CRBzOutput = ['z:  ',num2str(CRBinNanometer(3)), ' nm']; 
                end


                app.CRBOutputField.Text = [CRBxOutput newline CRByOutput newline CRBzOutput ];
            end

            app.CalculatingLamp.Color = "g";
        end

        function par = getParameters(app)
            % Fixed
            par.astigmatism = 0; % astigmatism coefficient (in units of wavelength; 0.11 corresponds to 75nm RMS wavefront error)
            
            % Variable
            par.wavelength = Length(app.EmissionWavelengthSpinner.Value,'nm');
            par.objectiveNA = app.ObjectiveNaSpinner.Value;
            par.refractiveIndices = [app.RefractiveIndexSpecimenSpinner.Value,...
                app.RefractiveIndexIntermediateLayerSpinner.Value,...
                app.RefractiveIndexImmersionMediumSpinner.Value];
            par.pixelSize = Length(app.PixelsizeObjectSpaceSpinner.Value,'nm'); % pixelsize in nm
            par.nPixels = ceil(1000.*app.ROIsidelengthEditField.Value/(par.pixelSize.inNanometer)); % show window approximately 2.5µm wide
            par.nPhotons = app.PhotonnumberSpinner.Value; % signal photon count
            par.backgroundNoise = (app.BackgroundnoisestdSpinner.Value)^2; % background mean photon count
            par.heightIntermediateLayer = Length(app.IntermediateLayerThicknessSpinner.Value, 'mu');
            
            if app.ShowPsf3DCheckBox.Value
                zStep = app.StepsizeEditField.Value;
                nSteps = app.NumberstepsEditField.Value - 1;
                par.defocus = Length(app.DefocusSlider.Value-nSteps/2*zStep:zStep:app.DefocusSlider.Value+nSteps/2*zStep,'nm'); % defocus in nm
            else
                par.defocus = Length(app.DefocusSlider.Value,'nm'); % defocus in nm
            end

            par.phaseMask = app.readParametersPhaseMask();
            par.transmission = app.readParametersTransmissionMask();
            [par.zernikeNollIndices, par.zernikeCoefficients] = app.readParametersZernikeAberrations();
            
            switch app.PhotonshotnoiseSwitch.Value
                case 'Off'
                    par.shotNoise = 0;
                case 'On'
                    par.shotNoise = 1;
            end
            switch app.ReducedexcitationSwitch.Value
                case 'Off'
                    par.reducedExcitation = 0;
                case 'On'
                    par.reducedExcitation = 1;
            end
        end

        function parPhaseMask = readParametersPhaseMask(app)
            phaseMask = app.BFPmanipulationDropDown.Value; % shape of BFP manipulation
            rotationAngle = app.RotatephasemaskSlider.Value * 2 *pi;
            innerRadius = app.InnerringradiusSlider.Value; % radius of hole in middle
            switch phaseMask
                case 'none'
                    parPhaseMask = @(n) EmptyPhaseMask(n);
                case 'Astigmatism'
                    parPhaseMask = @(n) Astigmatism(n,0.11).cutInnerRing(innerRadius).rotate(rotationAngle);
                case 'Vortex'
                    parPhaseMask = @(n) Vortex(n).cutInnerRing(innerRadius).rotate(rotationAngle);
                case 'Sector'
                    sector = deg2rad(app.SectorSlider.Value);
                    parPhaseMask = @(n) Sector(n,[0,sector]).cutInnerRing(innerRadius).rotate(rotationAngle);
                case 'OpposingSectors'
                    sector = deg2rad(app.SectorSlider.Value);
                    parPhaseMask = @(n) OpposingSectors(n,[0,sector]).cutInnerRing(innerRadius).rotate(rotationAngle);
                case 'Pyramid'
                    nFacets = app.NumberfacetsSpinner.Value;
                    maxShift = app.MaxShiftSlider.Value;
                    parPhaseMask = @(n) PyramidPhaseMask(n,nFacets,maxShift).cutInnerRing(innerRadius).rotate(rotationAngle);    
                case 'DoubleHelix'
                    parPhaseMask = @(n) DoubleHelix(n).cutInnerRing(innerRadius).rotate(rotationAngle);
                case 'Custom'
                    if isempty(app.PhaseMaskFilepathLabel.Text)
                        parPhaseMask = @(n) EmptyPhaseMask(n);
                    else
                        % Interpolate to get correct size for discretization of BFP
                        sz = size(app.phaseMask.mask);
                        xg = 1:sz(1);
                        yg = 1:sz(2);
                        F = griddedInterpolant({xg,yg},double(app.phaseMask.mask));
    
                        xq = linspace(1,sz(1),app.nDiscretizationBFP);
                        yq = linspace(1,sz(2),app.nDiscretizationBFP);
                        phaseMaskInterpolated = F({xq,yq});
                        parPhaseMask = PhaseMask(phaseMaskInterpolated);
                    end
            end
        end

        function parTransmissionMask = readParametersTransmissionMask(app)
            switch app.TransmissionmaskDropDown.Value
                case 'none'
                    parTransmissionMask = @(n) EmptyTransmissionMask(n);
                case 'Custom'
                    if isempty(app.TransmissionMaskFilepathLabel.Text)
                        parTransmissionMask = @(n) EmptyTransmissionMask(n);
                    else
                        % Interpolate to get correct size for discretization of BFP
                        sz = size(app.transmissionMask.mask);
                        xg = 1:sz(1);
                        yg = 1:sz(2);
                        F = griddedInterpolant({xg,yg},double(app.transmissionMask.mask));
    
                        xq = linspace(1,sz(1),app.nDiscretizationBFP);
                        yq = linspace(1,sz(2),app.nDiscretizationBFP);
                        transmissionMaskInterpolated = F({xq,yq});
                        parTransmissionMask = Transmission(transmissionMaskInterpolated);
                    end
            end
        end

        function [zernikeIndices, zernikeWeights] = readParametersZernikeAberrations(app)
            switch app.SpecifyZernikeaberrationsButtonGroup.SelectedObject.Text
                case 'Select common aberrations'
                    % Common aberrations specified
                    checkBoxes = {'VerticalTilt', 'HorizontalTilt', 'Defocus', 'PrimarySpherical', ...
                        'ObliqueAstigmatism', 'VerticalAstigmatism', 'VerticalComa', 'HorizontalComa'};
                    zernikeNollIndexing = [3, 2, 4, 11, 5, 6, 7, 8];
                    N = numel(checkBoxes);
                    
                    % Loop through all checkboxes
                    zernikeIndices = [];
                    zernikeWeights = [];
                    for i = 1:N
                        if app.([checkBoxes{i},'CheckBox']).Value
                            zernikeIndices = [zernikeIndices, zernikeNollIndexing(i)];
                            zernikeWeights = [zernikeWeights, app.(['Weight', checkBoxes{i}, 'Spinner']).Value];
                        end
                    end
                case 'Specify Zernike indices (Noll) and coefficients via input vector'
                    % Read from input field
                    if ~isempty(app.IndicesWeightsEditField.Value)
                        indicesText = app.IndicesWeightsEditField.Value;
                        indicesWeights = str2double(split(indicesText,{':',';'}));
                        indicesWeights = reshape(indicesWeights,2,[])';
                        zernikeIndices = indicesWeights(:,1);
                        zernikeWeights = indicesWeights(:,2);
                    else
                        zernikeIndices = [];
                        zernikeWeights = [];
                    end
                otherwise
                    error('Undefined option for SpecifyZernikeaberrationsButtonGroup.')
            end
            assert(length(zernikeWeights)==length(zernikeIndices))
            app.zernikeIndices = zernikeIndices;
            app.zernikeWeights = zernikeWeights;
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.originalPath = path; % save path at startup
            
            par.objectiveNA = 0.7;
            par.nPhotons = 50000;
            par.dipole = Dipole(0,0);
            par.position = Length([0 0 0],'nm');
            par.backgroundNoise = 0;
            par.defocus = Length(0,'nm');
            par.astigmatism = 0;
            par.reducedExcitation = 0;
            par.pixelSize = Length(100,'nm');
            par.nPixels = ceil(2500/(par.pixelSize.inNanometer)); % show window approximately 2.5µm wide 
            par.phaseMask = @(n) EmptyPhaseMask(n);
            
            % Create Zernike instance
            parMask.nGrid = 129;
            pupilMask = Mask(parMask);
            ZernikeInstance = ZernikePolynomials.getInstance(pupilMask);

            % Initialize phase mask
            app.phaseMask = EmptyPhaseMask(129);

            % Simulate PSF
            psf = PSF(par);
            psfImage = psf.image;

            % Open the PSF plot window and pass PSF image
            app.PlotPSF = WindowPlotPSF(app);
            app.PlotPSF.initializePlot(psfImage);
        end

        % Selection changed function: DipolerotationButtonGroup
        function DipolerotationButtonGroupSelectionChanged(app, event)
            switch app.DipolerotationButtonGroup.SelectedObject.Text
                case 'fixed'
                    if app.SwitchMultipleFluorophores
                        for k = 1:app.Fluorophores.getNumberFluorophores
                            app.Fluorophores.updateDipoleRotation('fixed',k);
                        end
                    else % single fluorophore case 
                        app.AzimuthalAngleSliderLabel.Visible = "on";
                        app.AzimuthalAngleSlider.Visible = "on";
                        app.phi.Visible = "on";
                        app.InclinationAngleSliderLabel.Visible = "on";
                        app.InclinationAngleSlider.Visible = "on";
                        app.theta.Visible = "on";
                        app.ReducedexcitationSwitch.Visible = "on";
                        app.ReducedexcitationSwitchLabel.Visible = "on";
                        app.ConfiguremultiplefluorophoresButton.Visible = "on";
                    end
                case 'freely rotating'
                    if app.SwitchMultipleFluorophores
                        for k = 1:app.Fluorophores.getNumberFluorophores
                            app.Fluorophores.updateDipoleRotation('freely rotating',k);
                        end
                    else
                        app.AzimuthalAngleSliderLabel.Visible = "off";
                        app.AzimuthalAngleSlider.Visible = "off";
                        app.phi.Visible = "off";
                        app.InclinationAngleSliderLabel.Visible = "off";
                        app.InclinationAngleSlider.Visible = "off";
                        app.theta.Visible = "off";
                        app.ReducedexcitationSwitch.Visible = "off";
                        app.ReducedexcitationSwitchLabel.Visible = "off";
                        app.ConfiguremultiplefluorophoresButton.Visible = "on";
                    end
                otherwise
                    error('Invalid input value for dipole rotation!')
            end
            
            app.DipolerotationButtonGroup.SelectedObject = app.DipolerotationButtonGroup.SelectedObject;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: PhotonnumberSpinner
        function PhotonnumberSpinnerValueChanging(app, event)
            if ischar(event.Value)
                value = str2double(event.Value);
            else
                value = event.Value;
            end
            app.PhotonnumberSpinner.Value = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: PhotonnumberSpinner
        function PhotonnumberSpinnerValueChanged(app, event)
            app.PhotonnumberSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: ReducedexcitationSwitch
        function ReducedexcitationSwitchValueChanged(app, event)
            app.ReducedexcitationSwitch.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: PhotonshotnoiseSwitch
        function PhotonshotnoiseSwitchValueChanged(app, event)
            app.PhotonshotnoiseSwitch.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: InclinationAngleSlider
        function InclinationAngleSliderValueChanging(app, event)
            app.InclinationAngleSlider.Value = event.Value;
            app.theta.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: InclinationAngleSlider
        function InclinationAngleSliderValueChanged(app, event)
            app.InclinationAngleSlider.Value = event.Value;
            app.theta.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: theta
        function thetaValueChanged(app, event)
            app.theta.Value = event.Value;
            app.InclinationAngleSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: AzimuthalAngleSlider
        function AzimuthalAngleSliderValueChanging(app, event)
            app.AzimuthalAngleSlider.Value = event.Value;
            app.phi.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: AzimuthalAngleSlider
        function AzimuthalAngleSliderValueChanged(app, event)
            app.AzimuthalAngleSlider.Value = event.Value;
            app.phi.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: phi
        function phiValueChanged(app, event)
            app.phi.Value = event.Value;
            app.AzimuthalAngleSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: zpositionSpinner
        function zpositionSpinnerValueChanging(app, event)
            app.zpositionSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: zpositionSpinner
        function zpositionSpinnerValueChanged(app, event)
            app.zpositionSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Button pushed function: ConfiguremultiplefluorophoresButton
        function ConfiguremultiplefluorophoresButtonPushed(app, event)
            % Multiple fluorophores currently not compatible with 3D plot
            if app.ShowPsf3DCheckBox.Value == true
                app.ShowPsf3DCheckBox.Value = false;
                app.ShowPsf3DCheckBoxValueChanged()
                app.ShowPsf3DCheckBox.Enable = "off";
            end
            
            app.SwitchMultipleFluorophores = 1;

            app.AzimuthalAngleSliderLabel.Visible = "off";
            app.AzimuthalAngleSlider.Visible = "off";
            app.phi.Visible = "off";
            app.InclinationAngleSliderLabel.Visible = "off";
            app.InclinationAngleSlider.Visible = "off";
            app.theta.Visible = "off";

            app.xpositionSpinnerLabel.Visible = "off";
            app.xpositionSpinner.Visible = "off";
            app.ypositionSpinnerLabel.Visible = "off";
            app.ypositionSpinner.Visible = "off";
            app.zpositionSpinnerLabel.Visible = "off";
            app.zpositionSpinner.Visible = "off";

            app.CalculateCramrRaoBoundCheckBox.Visible = 'off'; 
            app.CRBOutputField.Visible = 'off';
            app.CramrRaoBoundLabel.Visible = 'off';
          
            app.Fluorophores = WindowFluorophores(app);
            set(app.ConfiguremultiplefluorophoresButton, 'Enable', 'off')
        end

        % Value changing function: ObjectiveNaSpinner
        function ObjectiveNaSpinnerValueChanging(app, event)
            app.ObjectiveNaSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: ObjectiveNaSpinner
        function ObjectiveNaSpinnerValueChanged(app, event)
            app.ObjectiveNaSpinner.Value = app.ObjectiveNaSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: DefocusSlider
        function DefocusSliderValueChanging(app, event)
            app.DefocusSlider.Value = event.Value;
            app.defocus.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: DefocusSlider
        function DefocusSliderValueChanged(app, event)
            app.DefocusSlider.Value = event.Value;
            app.defocus.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: defocus
        function defocusValueChanged(app, event)
            app.defocus.Value = event.Value;
            app.DefocusSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Button pushed function: FocusButton
        function FocusButtonPushed(app, event)
            app.DefocusSlider.Value = 0;
            app.defocus.Value = 0;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: BackgroundnoisestdSpinner
        function BackgroundnoisestdSpinnerValueChanging(app, event)
            if ischar(event.Value)
                value = str2double(event.Value);
            else
                value = event.Value;
            end
            app.BackgroundnoisestdSpinner.Value = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: BackgroundnoisestdSpinner
        function BackgroundnoisestdSpinnerValueChanged(app, event)
            app.BackgroundnoisestdSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: BFPmanipulationDropDown
        function BFPmanipulationDropDownValueChanged(app, event)
            value = event.Value;
            app.BFPmanipulationDropDown.Value = value;
            switch value
                case {'none','Custom'}
                    app.PhaseMaskOptionsLabel.Visible = "off";
                    app.InnerringradiusSliderLabel.Visible = "off";
                    app.InnerringradiusSlider.Visible = "off";
                    app.RotatephasemaskSliderLabel.Visible = "off";
                    app.RotatephasemaskSlider.Visible = "off";
                otherwise
                    app.PhaseMaskOptionsLabel.Visible = "on";
                    app.InnerringradiusSliderLabel.Visible = "on";
                    app.InnerringradiusSlider.Visible = "on";
                    app.RotatephasemaskSliderLabel.Visible = "on";
                    app.RotatephasemaskSlider.Visible = "on";
            end
            switch value
                case 'Sector'
                    app.SectorSliderLabel.Visible = "on";
                    app.SectorSlider.Visible = "on";
                    app.SectorSlider.Limits = [0,360];
                    app.SectorSlider.Value = 180;
                    app.NumberfacetsSpinnerLabel.Visible = "off";
                    app.NumberfacetsSpinner.Visible = "off";
                    app.MaxShiftSliderLabel.Visible = "off";
                    app.MaxShiftSlider.Visible = "off";
                case 'OpposingSectors'
                    app.SectorSliderLabel.Visible = "on";
                    app.SectorSlider.Visible = "on";
                    app.SectorSlider.Limits = [0,180];
                    app.SectorSlider.Value = 90;
                    app.NumberfacetsSpinnerLabel.Visible = "off";
                    app.NumberfacetsSpinner.Visible = "off";
                    app.MaxShiftSliderLabel.Visible = "off";
                    app.MaxShiftSlider.Visible = "off";
                case 'Pyramid'
                    app.SectorSliderLabel.Visible = "off";
                    app.SectorSlider.Visible = "off";
                    app.NumberfacetsSpinnerLabel.Visible = "on";
                    app.NumberfacetsSpinner.Visible = "on";
                    app.MaxShiftSliderLabel.Visible = "on";
                    app.MaxShiftSlider.Visible = "on";
                otherwise
                    app.SectorSliderLabel.Visible = "off";
                    app.SectorSlider.Visible = "off";
                    app.NumberfacetsSpinnerLabel.Visible = "off";
                    app.NumberfacetsSpinner.Visible = "off";
                    app.MaxShiftSliderLabel.Visible = "off";
                    app.MaxShiftSlider.Visible = "off";
            end
            switch value
                case 'Custom'
                    app.LoadCustomPhaseMaskButton.Visible = "on";
                    app.PhaseMaskFilepathLabel.Visible = "on";
                otherwise
                    app.LoadCustomPhaseMaskButton.Visible = "off";
                    app.PhaseMaskFilepathLabel.Visible = "off";
            end
            simulateAndDisplayPSF(app);
        end

        % Value changing function: InnerringradiusSlider
        function InnerringradiusSliderValueChanging(app, event)
            app.InnerringradiusSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: InnerringradiusSlider
        function InnerringradiusSliderValueChanged(app, event)
            app.InnerringradiusSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: RotatephasemaskSlider
        function RotatephasemaskSliderValueChanging(app, event)
            app.RotatephasemaskSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: RotatephasemaskSlider
        function RotatephasemaskSliderValueChanged(app, event)
            app.RotatephasemaskSlider.Value = app.RotatephasemaskSlider.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: SectorSlider
        function SectorSliderValueChanging(app, event)
            app.SectorSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: SectorSlider
        function SectorSliderValueChanged(app, event)
            app.SectorSlider.Value = app.SectorSlider.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: NumberfacetsSpinner
        function NumberfacetsSpinnerValueChanging(app, event)
            app.NumberfacetsSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: NumberfacetsSpinner
        function NumberfacetsSpinnerValueChanged(app, event)
            app.NumberfacetsSpinner.Value = app.NumberfacetsSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: MaxShiftSlider
        function MaxShiftSliderValueChanging(app, event)
            app.MaxShiftSlider.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: MaxShiftSlider
        function MaxShiftSliderValueChanged(app, event)
            app.MaxShiftSlider.Value = app.MaxShiftSlider.Value;
            simulateAndDisplayPSF(app);
        end

        % Close request function: PSFsimulationUIFigure
        function PSFsimulationUIFigureCloseRequest(app, event)
            delete(app.PlotPSF)
            delete(app.PlotPSFThreeDim)
            delete(app.PlotPhaseMask)
            delete(app.PlotTransmissionMask)
            delete(app.PlotZernikeAberrations)
            delete(app.Fluorophores)
            delete(app.PlotPolarizedEmission)
            delete(app.FitZernikeApp)
            delete(app)
        end

        % Selection changed function: SetcontrastButtonGroup
        function SetcontrastButtonGroupSelectionChanged(app, event)
            app.SetcontrastButtonGroup.SelectedObject = event.NewValue;

            if app.ShowPSFCheckBox.Value
                app.PlotPSF.setContrast( app.SetcontrastButtonGroup.SelectedObject.Text );
            end
            if app.SwitchMultipleFluorophores == 1
                app.Fluorophores.setContrastForAll( app.SetcontrastButtonGroup.SelectedObject.Text );
            end
        end

        % Value changed function: ZernikeAberrationsShowplotCheckBox
        function ZernikeAberrationsShowplotCheckBoxValueChanged(app, event)
            value = app.ZernikeAberrationsShowplotCheckBox.Value;
            if value
                app.PlotZernikeAberrations = WindowZernikeAberrations(app);
                app.PlotZernikeAberrations.initializePlot();
            else
                delete(app.PlotZernikeAberrations);
            end
        end

        % Value changed function: PhaseMaskShowplotCheckBox
        function PhaseMaskShowplotCheckBoxValueChanged(app, event)
            value = app.PhaseMaskShowplotCheckBox.Value;
            if value
                app.PlotPhaseMask = WindowPhaseMask(app);
                app.PlotPhaseMask.initializePlot();
            else
                delete(app.PlotPhaseMask);
            end
        end

        % Value changed function: ShowPSFCheckBox
        function ShowPSFCheckBoxValueChanged(app, event)
            value = app.ShowPSFCheckBox.Value;
            if value
                app.CalculatingLamp.Enable = 'on';
                app.ShowPsf2DCheckBox.Value = true;
                app.ShowPsf2DCheckBox.Enable = "on";
                app.ShowPsf3DCheckBox.Enable = "on";
                app.PlotPSF = WindowPlotPSF(app);
                app.PlotPSF.initializePlot();
            else
                appPath = app.originalPath;
                delete(app.PlotPSF);
                path(appPath);
                app.CalculatingLamp.Enable = 'off';
                app.ShowPsf2DCheckBox.Value = false;
                app.ShowPsf3DCheckBox.Value = false;
                app.ShowPsf2DCheckBox.Enable = "off";
                app.ShowPsf3DCheckBox.Enable = "off";
            end
        end

        % Value changed function: ShowPsf2DCheckBox
        function ShowPsf2DCheckBoxValueChanged(app, event)
            value = app.ShowPsf2DCheckBox.Value;
            if value
                app.PlotPSF = WindowPlotPSF(app);
                app.PlotPSF.initializePlot();
            else
                appPath = app.originalPath;
                delete(app.PlotPSF);
                path(appPath);
                if app.ShowPsf3DCheckBox.Value == false
                    app.ShowPSFCheckBox.Value = false;
                    app.ShowPsf2DCheckBox.Enable = "off";
                    app.ShowPsf3DCheckBox.Enable = "off";
                end
            end
        end

        % Value changed function: ShowPsf3DCheckBox
        function ShowPsf3DCheckBoxValueChanged(app, event)
            value = app.ShowPsf3DCheckBox.Value;
            if value
                app.PlotPSFThreeDim = WindowPlotPSFThreeDim(app);
                app.PlotPSFThreeDim.initializePlot();
                app.StepsizeEditFieldLabel.Visible = "on";
                app.StepsizeEditFieldLabel.Enable = "on";
                app.StepsizeEditField.Visible = "on";
                app.StepsizeEditField.Enable = "on";
                app.NumberstepsEditFieldLabel.Visible = "on";
                app.NumberstepsEditFieldLabel.Enable = "on";
                app.NumberstepsEditField.Visible = "on";
                app.NumberstepsEditField.Enable = "on";
                simulateAndDisplayPSF(app);
            else
                delete(app.PlotPSFThreeDim);
                app.StepsizeEditFieldLabel.Visible = "off";
                app.StepsizeEditFieldLabel.Enable = "off";
                app.StepsizeEditField.Visible = "off";
                app.StepsizeEditField.Enable = "off";
                app.NumberstepsEditFieldLabel.Visible = "off";
                app.NumberstepsEditFieldLabel.Enable = "off";
                app.NumberstepsEditField.Visible = "off";
                app.NumberstepsEditField.Enable = "off";
                if app.ShowPsf2DCheckBox.Value == false
                    app.ShowPSFCheckBox.Value = false;
                    app.ShowPsf2DCheckBox.Enable = "off";
                    app.ShowPsf3DCheckBox.Enable = "off";
                end
            end
        end

        % Value changed function: PolarizedemissionchannelsCheckBox
        function PolarizedemissionchannelsCheckBoxValueChanged(app, event)
            app.PolarizedemissionchannelsCheckBox.Value = event.Value;
            if app.PolarizedemissionchannelsCheckBox.Value
                % Open the Plot window and pass psf object
                app.PlotPolarizedEmission = WindowPolarizedEmission(app);
                app.PlotPolarizedEmission.initializePlot();
                switch app.SwitchMultipleFluorophores
                    case 0 % 'Single'
                        psf = simulateAndDisplayPSF(app);
                        app.PlotPolarizedEmission.updatePlot(psf);
                    case 1 % 'Multiple'
                        app.Fluorophores.updatePolarizedEmissionChannels();
                end
            else
                delete(app.PlotPolarizedEmission);
            end
        end

        % Value changed function: VerticalTiltCheckBox
        function VerticalTiltCheckBoxValueChanged(app, event)
            value = app.VerticalTiltCheckBox.Value;
            app.WeightVerticalTiltSpinner.Enable = value;
            app.WeightVerticalTiltSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: HorizontalTiltCheckBox
        function HorizontalTiltCheckBoxValueChanged(app, event)
            value = app.HorizontalTiltCheckBox.Value;
            app.WeightHorizontalTiltSpinner.Enable = value;
            app.WeightHorizontalTiltSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: DefocusCheckBox
        function DefocusCheckBoxValueChanged(app, event)
            value = app.DefocusCheckBox.Value;
            app.WeightDefocusSpinner.Enable = value;
            app.WeightDefocusSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: PrimarySphericalCheckBox
        function PrimarySphericalCheckBoxValueChanged(app, event)
            value = app.PrimarySphericalCheckBox.Value;
            app.WeightPrimarySphericalSpinner.Enable = value;
            app.WeightPrimarySphericalSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: ObliqueAstigmatismCheckBox
        function ObliqueAstigmatismCheckBoxValueChanged(app, event)
            value = app.ObliqueAstigmatismCheckBox.Value;
            app.WeightObliqueAstigmatismSpinner.Enable = value;
            app.WeightObliqueAstigmatismSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: VerticalAstigmatismCheckBox
        function VerticalAstigmatismCheckBoxValueChanged(app, event)
            value = app.VerticalAstigmatismCheckBox.Value;
            app.WeightVerticalAstigmatismSpinner.Enable = value;
            app.WeightVerticalAstigmatismSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: VerticalComaCheckBox
        function VerticalComaCheckBoxValueChanged(app, event)
            value = app.VerticalComaCheckBox.Value;
            app.WeightVerticalComaSpinner.Enable = value;
            app.WeightVerticalComaSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: HorizontalComaCheckBox
        function HorizontalComaCheckBoxValueChanged(app, event)
            value = app.HorizontalComaCheckBox.Value;
            app.WeightHorizontalComaSpinner.Enable = value;
            app.WeightHorizontalComaSpinnerLabel.Enable = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightVerticalTiltSpinner
        function WeightVerticalTiltSpinnerValueChanged(app, event)
            app.WeightVerticalTiltSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightHorizontalTiltSpinner
        function WeightHorizontalTiltSpinnerValueChanged(app, event)
            app.WeightHorizontalTiltSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightDefocusSpinner
        function WeightDefocusSpinnerValueChanged(app, event)
            app.WeightDefocusSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightPrimarySphericalSpinner
        function WeightPrimarySphericalSpinnerValueChanged(app, event)
            app.WeightPrimarySphericalSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightObliqueAstigmatismSpinner
        function WeightObliqueAstigmatismSpinnerValueChanged(app, event)
            app.WeightObliqueAstigmatismSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightVerticalAstigmatismSpinner
        function WeightVerticalAstigmatismSpinnerValueChanged(app, event)
            app.WeightVerticalAstigmatismSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightVerticalComaSpinner
        function WeightVerticalComaSpinnerValueChanged(app, event)
            app.WeightVerticalComaSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: WeightHorizontalComaSpinner
        function WeightHorizontalComaSpinnerValueChanged(app, event)
            app.WeightHorizontalComaSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Button pushed function: FitZernikeButton
        function FitZernikeButtonPushed(app, event)
            app.FitZernikeApp = aberration_measurement(app);
        end

        % Button pushed function: LoadCustomPhaseMaskButton
        function LoadCustomPhaseMaskButtonPushed(app, event)
            [filename, pathname] = uigetfile('*.mat');
            try
                mask = loadSingleMatrix(fullfile(pathname, filename));
                app.phaseMask = PhaseMask(mask);
                if ~isempty(app.PlotPhaseMask)
                    app.PlotPhaseMask.updatePlot();
                end
                simulateAndDisplayPSF(app);
                app.PhaseMaskFilepathLabel.Text = ['Loaded file: ', filename];
            catch
                warndlg('Error loading file.','Warning');
            end
        end

        % Button pushed function: LoadZernikeIndicesWeightsButton
        function LoadZernikeIndicesWeightsButtonPushed(app, event)
            [filename, pathname] = uigetfile({'*.mat';'*.csv'},'Select file for Zernike indices and coefficients');
            [app.zernikeIndices, app.zernikeWeights] = loadZernikeIndicesAndWeights(fullfile(pathname, filename));
            
            indicesCell = cellstr(num2str(app.zernikeIndices(:)));
            weightsCell = cellstr(num2str(app.zernikeWeights(:)));
            indicesWeightsCell = strcat(indicesCell, ':', weightsCell);
            
            app.IndicesWeightsEditField.Value = strjoin(indicesWeightsCell, ';');
            if ~isempty(app.IndicesWeightsEditField.Value)
                app.SaveZernikeIndicesWeightsButton.Enable = 'on';
            else
                app.SaveZernikeIndicesWeightsButton.Enable = 'off';
            end
            if ~isempty(app.PlotZernikeAberrations)
                app.PlotZernikeAberrations.updatePlot();
            end
            simulateAndDisplayPSF(app);
        end

        % Selection changed function: SpecifyZernikeaberrationsButtonGroup
        function SpecifyZernikeaberrationsButtonGroupSelectionChanged(app, event)
            selectedButton = app.SpecifyZernikeaberrationsButtonGroup.SelectedObject;
            handleArrayCommonAberrationsCheckBoxes = [app.VerticalTiltCheckBox, app.HorizontalTiltCheckBox,...
                app.DefocusCheckBox, app.PrimarySphericalCheckBox,...
                app.ObliqueAstigmatismCheckBox, app.VerticalAstigmatismCheckBox,...
                app.VerticalComaCheckBox, app.HorizontalComaCheckBox];
            handleArrayCommonAberrationsWeights = [app.WeightHorizontalComaSpinner, app.WeightVerticalComaSpinner,...
                app.WeightVerticalAstigmatismSpinner, app.WeightObliqueAstigmatismSpinner,...
                app.WeightPrimarySphericalSpinner, app.WeightDefocusSpinner,...
                app.WeightHorizontalTiltSpinner, app.WeightVerticalTiltSpinner,...
                app.WeightHorizontalComaSpinnerLabel, app.WeightVerticalComaSpinnerLabel,...
                app.WeightVerticalAstigmatismSpinnerLabel, app.WeightObliqueAstigmatismSpinnerLabel,...
                app.WeightPrimarySphericalSpinnerLabel, app.WeightDefocusSpinnerLabel,...
                app.WeightHorizontalTiltSpinnerLabel, app.WeightVerticalTiltSpinnerLabel];
            handleArrayInputVector = [app.LoadZernikeIndicesWeightsButton, app.IndicesWeightsEditField];
            switch selectedButton.Text
                case 'Select common aberrations'
                    set(handleArrayInputVector, 'Enable', 'off');
                    set(handleArrayCommonAberrationsCheckBoxes, 'Enable', 'on');
                    set(handleArrayCommonAberrationsCheckBoxes, 'Value', false);
                case 'Specify Zernike indices (Noll) and coefficients via input vector'
                    set(handleArrayCommonAberrationsCheckBoxes, 'Enable', 'off');
                    set(handleArrayCommonAberrationsWeights, 'Enable', 'off');
                    set(handleArrayCommonAberrationsCheckBoxes, 'Value', false);
                    set(handleArrayInputVector, 'Enable', 'on');
                otherwise
                    error('Undefined option for SpecifyZernikeaberrationsButtonGroup.')
            end
            if ~isempty(app.PlotZernikeAberrations)
                app.PlotZernikeAberrations.updatePlot();
            end
            simulateAndDisplayPSF(app);
        end

        % Value changed function: IndicesWeightsEditField
        function IndicesWeightsEditFieldValueChanged(app, event)
            indicesText = app.IndicesWeightsEditField.Value;
            % Assert that string is formatted correctly
            pattern = '^\d+:\d*\.?\d+(;\d+:\d*\.?\d+)*$';
            if isempty(regexp(indicesText, pattern, 'once'))
                uialert(app.PSFsimulationUIFigure,'Invalid input! Specify Zernike Noll indices and coefficients in the correct format.','Invalid input');
            else
                indicesWeights = str2double(split(indicesText,{':',';'}));
                indicesWeights = reshape(indicesWeights,2,[])';
                indices = indicesWeights(:,1);
                if any(indices==0) % Check for 0
                    uialert(app.PSFsimulationUIFigure,'Invalid input! Index 0 is not a valid Zernike Noll index.','Invalid input');
                elseif numel(indices)~=numel(unique(indices)) % Check for duplicated indices
                    uialert(app.PSFsimulationUIFigure,'Invalid input! Duplicated Zernike indices.','Invalid input');
                else
                    app.zernikeIndices = indicesWeights(:,1);
                    app.zernikeWeights = indicesWeights(:,2);
                    if ~isempty(app.PlotZernikeAberrations)
                        app.PlotZernikeAberrations.updatePlot();
                    end
                    simulateAndDisplayPSF(app);
                end
            end
            if ~isempty(app.IndicesWeightsEditField.Value)
                app.SaveZernikeIndicesWeightsButton.Enable = 'on';
            else
                app.SaveZernikeIndicesWeightsButton.Enable = 'off';
            end
        end

        % Value changed function: RefractiveIndexImmersionMediumSpinner
        function RefractiveIndexImmersionMediumSpinnerValueChanged(app, event)
            app.RefractiveIndexImmersionMediumSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: RefractiveIndexIntermediateLayerSpinner
        function RefractiveIndexIntermediateLayerSpinnerValueChanged(app, event)
            app.RefractiveIndexIntermediateLayerSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: RefractiveIndexSpecimenSpinner
        function RefractiveIndexSpecimenSpinnerValueChanged(app, event)
            app.RefractiveIndexSpecimenSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: RefractiveIndexImmersionMediumSpinner
        function RefractiveIndexImmersionMediumSpinnerValueChanging(app, event)
            app.RefractiveIndexImmersionMediumSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changing function: RefractiveIndexIntermediateLayerSpinner
        function RefractiveIndexIntermediateLayerSpinnerValueChanging(app, event)
            app.RefractiveIndexIntermediateLayerSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Button pushed function: SaveZernikeIndicesWeightsButton
        function SaveZernikeIndicesWeightsButtonPushed(app, event)
            zernike = [app.zernikeIndices, app.zernikeWeights];
            
            [filename, pathname] = uiputfile({'*.mat';'*.csv'}, 'Save Zernike indices and coefficients');
            fullFilename = fullfile(pathname, filename);
            [~,~,extension] = fileparts(fullFilename);
            
            if filename ~= 0 % if user did not cancel the file dialog
                switch extension
                    case '.mat'
                        save(fullFilename,zernike,'-mat')
                    case '.csv'
                        % write the matrix to the CSV file
                        writematrix(zernike, fullFilename);
                    otherwise
                        error('File extension not implemented!')
                end
            end
        end

        % Value changed function: EmissionWavelengthSpinner
        function EmissionWavelengthSpinnerValueChanged(app, event)
            app.EmissionWavelengthSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: ObjectiveFocalLengthSpinner
        function ObjectiveFocalLengthSpinnerValueChanged(app, event)
            app.ObjectiveFocalLengthSpinner.Value = event.Value;
            app.MagnificationSpinner.Value = app.TubeLensFocalLengthSpinner.Value / app.ObjectiveFocalLengthSpinner.Value;
            app.PixelsizeObjectSpaceSpinner.Value = app.PixelsizePhysicalSpinner.Value*1e3 / app.MagnificationSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: TubeLensFocalLengthSpinner
        function TubeLensFocalLengthSpinnerValueChanged(app, event)
            app.TubeLensFocalLengthSpinner.Value = event.Value;
            app.MagnificationSpinner.Value = app.TubeLensFocalLengthSpinner.Value / app.ObjectiveFocalLengthSpinner.Value;
            app.PixelsizeObjectSpaceSpinner.Value = app.PixelsizePhysicalSpinner.Value*1e3 / app.MagnificationSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: MagnificationSpinner
        function MagnificationSpinnerValueChanged(app, event)
            app.MagnificationSpinner.Value = event.Value;
            app.ObjectiveFocalLengthSpinner.Value = app.TubeLensFocalLengthSpinner.Value / app.MagnificationSpinner.Value;
            app.PixelsizeObjectSpaceSpinner.Value = app.PixelsizePhysicalSpinner.Value*1e3 / app.MagnificationSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: PixelsizeObjectSpaceSpinner
        function PixelsizeObjectSpaceSpinnerValueChanged(app, event)
            app.PixelsizeObjectSpaceSpinner.Value = event.Value;
            app.PixelsizePhysicalSpinner.Value = app.PixelsizeObjectSpaceSpinner.Value*1e-3 * app.MagnificationSpinner.Value;
            
            pixelSize = event.Value;
            nPixels = ceil(app.ROIsidelengthEditField.Value/pixelSize*1e3); % ROI size given in µm, pixelSize in nm
            roundedValue = nPixels*pixelSize*1e-3;
            if roundedValue > app.ROIsidelengthEditField.Limits(2)
                roundedValue = roundedValue - 2*pixelSize*1e-3;
            end
            app.ROIsidelengthEditField.Value = roundedValue;
            app.nPixels = nPixels;
            app.PixelsperlateralaxisEditField.Value = nPixels;

            simulateAndDisplayPSF(app);
        end

        % Value changing function: PixelsizePhysicalSpinner
        function PixelsizePhysicalSpinnerValueChanging(app, event)
            if ischar(event.Value)
                value = str2double(event.Value);
            else
                value = event.Value;
            end
            app.PixelsizePhysicalSpinner.Value = value;
            app.PixelsizeObjectSpaceSpinner.Value = app.PixelsizePhysicalSpinner.Value*1e3 / app.MagnificationSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: PixelsizePhysicalSpinner
        function PixelsizePhysicalSpinnerValueChanged(app, event)
            app.PixelsizePhysicalSpinner.Value = event.Value;
            app.PixelsizeObjectSpaceSpinner.Value = app.PixelsizePhysicalSpinner.Value*1e3 / app.MagnificationSpinner.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: AddintermediatelayerCheckBox
        function AddintermediatelayerCheckBoxValueChanged(app, event)
            app.AddintermediatelayerCheckBox.Value = event.Value;
            flagOnOff = app.AddintermediatelayerCheckBox.Value;
            app.RefractiveIndexIntermediateLayerSpinner.Visible = flagOnOff;
            app.RefractiveIndexIntermediateLayerLabel.Visible = flagOnOff;
            app.IntermediateLayerThicknessSpinner.Visible = flagOnOff;
            app.IntermediateLayerThicknessLabel.Visible = flagOnOff;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: IntermediateLayerThicknessSpinner
        function IntermediateLayerThicknessSpinnerValueChanged(app, event)
            app.IntermediateLayerThicknessSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: ColormapDropDown
        function ColormapDropDownValueChanged(app, event)
            value = event.Value;
            app.ColormapDropDown.Value = value;
            if app.ShowPSFCheckBox.Value && app.ShowPsf2DCheckBox.Value
                app.PlotPSF.setColormap(value);
            end
            if app.ShowPSFCheckBox.Value && app.ShowPsf3DCheckBox.Value
                app.PlotPSFThreeDim.setColormap(value);
            end
            if app.PolarizedemissionchannelsCheckBox.Value
                app.PlotPolarizedEmission.setColormap(value);
            end
            if app.SwitchMultipleFluorophores
                app.Fluorophores.setColormap(value)
            end
        end

        % Value changed function: ROIsidelengthEditField
        function ROIsidelengthEditFieldValueChanged(app, event)
            pixelSize = app.PixelsizeObjectSpaceSpinner.Value;

            nPixels = ceil(event.Value/pixelSize*1e3); % ROI size given in µm, pixelSize in nm
            roundedValue = nPixels*pixelSize*1e-3; 
            
            if roundedValue > app.ROIsidelengthEditField.Limits(2)
                roundedValue = roundedValue - 2*pixelSize*1e-3;
            end

            app.ROIsidelengthEditField.Value = roundedValue;
            app.nPixels = nPixels;
            app.PixelsperlateralaxisEditField.Value = nPixels;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: PixelsperlateralaxisEditField
        function PixelsperlateralaxisEditFieldValueChanged(app, event)
            nPixels = event.Value;
            app.PixelsperlateralaxisEditField.Value = nPixels;
            app.nPixels = nPixels;
            app.ROIsidelengthEditField.Value = double(app.nPixels) * app.PixelsizeObjectSpaceSpinner.Value * 1e-3;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: TransmissionmaskDropDown
        function TransmissionmaskDropDownValueChanged(app, event)
            value = event.Value;
            app.TransmissionmaskDropDown.Value = value;
            switch value
                case 'Custom'
                    app.LoadCustomTransmissionMaskButton.Visible = "on";
                    app.TransmissionMaskFilepathLabel.Visible = "on";
                otherwise
                    app.LoadCustomTransmissionMaskButton.Visible = "off";
                    app.TransmissionMaskFilepathLabel.Visible = "off";
            end
            if app.isTransmissionMaskLoaded
                simulateAndDisplayPSF(app);
            end
        end

        % Value changed function: TransmissionMaskShowplotCheckBox
        function TransmissionMaskShowplotCheckBoxValueChanged(app, event)
            value = app.TransmissionMaskShowplotCheckBox.Value;
            if value
                app.PlotTransmissionMask = WindowTransmissionMask(app);
                app.PlotTransmissionMask.initializePlot();
            else
                delete(app.PlotTransmissionMask);
            end
        end

        % Button pushed function: LoadCustomTransmissionMaskButton
        function LoadCustomTransmissionMaskButtonPushed(app, event)
            [filename, pathname] = uigetfile('*.mat');
            try
                mask = loadSingleMatrix(fullfile(pathname, filename));
                app.transmissionMask = Transmission(mask);
                if ~isempty(app.PlotTransmissionMask)
                    app.PlotTransmissionMask.updatePlot();
                end
                simulateAndDisplayPSF(app);
                app.TransmissionMaskFilepathLabel.Text = ['Loaded file: ', filename];
                app.isTransmissionMaskLoaded = true;
            catch
                warndlg('Error loading file.','Warning');
            end
        end

        % Value changed function: NumberstepsEditField
        function NumberstepsEditFieldValueChanged(app, event)
            value = event.Value;
            app.NumberstepsEditField.Value = value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: StepsizeEditField
        function StepsizeEditFieldValueChanged(app, event)
            app.StepsizeEditField.Value = event.Value;
            simulateAndDisplayPSF(app);
        end

        % Value changed function: CalculateCramrRaoBoundCheckBox
        function CalculateCramrRaoBoundCheckBoxValueChanged(app, event)
            app.CalculateCramrRaoBoundCheckBox.Value = event.Value;
            if event.Value
                simulateAndDisplayPSF(app);
            else
                app.CRBOutputField.Text = ''; 
            end
        end

        % Selection change function: TabGroup
        function TabGroupSelectionChanged(app, event)
            % calculate CRB if Options tab is selected
            if strcmp(app.TabGroup.SelectedTab.Title, 'Options')
                simulateAndDisplayPSF(app);
            end
        end

        % Value changed function: xpositionSpinner
        function xpositionSpinnerValueChanged(app, event)
            app.xpositionSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
            
        end

        % Value changing function: xpositionSpinner
        function xpositionSpinnerValueChanging(app, event)
            app.xpositionSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
            
        end

        % Value changed function: ypositionSpinner
        function ypositionSpinnerValueChanged(app, event)
            app.ypositionSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
            
        end

        % Value changing function: ypositionSpinner
        function ypositionSpinnerValueChanging(app, event)
            app.ypositionSpinner.Value = event.Value;
            simulateAndDisplayPSF(app);
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PSFsimulationUIFigure and hide until all components are created
            app.PSFsimulationUIFigure = uifigure('Visible', 'off');
            app.PSFsimulationUIFigure.Position = [100 120 633 469];
            app.PSFsimulationUIFigure.Name = 'PSF simulation';
            app.PSFsimulationUIFigure.CloseRequestFcn = createCallbackFcn(app, @PSFsimulationUIFigureCloseRequest, true);

            % Create SetparametersLabel
            app.SetparametersLabel = uilabel(app.PSFsimulationUIFigure);
            app.SetparametersLabel.FontSize = 14;
            app.SetparametersLabel.FontWeight = 'bold';
            app.SetparametersLabel.Position = [14 434 107 22];
            app.SetparametersLabel.Text = 'Set parameters';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.PSFsimulationUIFigure);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Position = [14 14 612 411];

            % Create FluorophoreTab
            app.FluorophoreTab = uitab(app.TabGroup);
            app.FluorophoreTab.Title = 'Fluorophore';

            % Create FluorophoreLabel
            app.FluorophoreLabel = uilabel(app.FluorophoreTab);
            app.FluorophoreLabel.FontSize = 13;
            app.FluorophoreLabel.FontWeight = 'bold';
            app.FluorophoreLabel.Position = [12 354 82 22];
            app.FluorophoreLabel.Text = 'Fluorophore';

            % Create AzimuthalAngleSliderLabel
            app.AzimuthalAngleSliderLabel = uilabel(app.FluorophoreTab);
            app.AzimuthalAngleSliderLabel.Position = [273 74 92 22];
            app.AzimuthalAngleSliderLabel.Text = 'Azimuthal Angle';

            % Create AzimuthalAngleSlider
            app.AzimuthalAngleSlider = uislider(app.FluorophoreTab);
            app.AzimuthalAngleSlider.Limits = [0 360];
            app.AzimuthalAngleSlider.ValueChangedFcn = createCallbackFcn(app, @AzimuthalAngleSliderValueChanged, true);
            app.AzimuthalAngleSlider.ValueChangingFcn = createCallbackFcn(app, @AzimuthalAngleSliderValueChanging, true);
            app.AzimuthalAngleSlider.BusyAction = 'cancel';
            app.AzimuthalAngleSlider.Position = [394 84 150 3];

            % Create InclinationAngleSliderLabel
            app.InclinationAngleSliderLabel = uilabel(app.FluorophoreTab);
            app.InclinationAngleSliderLabel.Position = [273 121 94 22];
            app.InclinationAngleSliderLabel.Text = 'Inclination Angle';

            % Create InclinationAngleSlider
            app.InclinationAngleSlider = uislider(app.FluorophoreTab);
            app.InclinationAngleSlider.Limits = [0 90];
            app.InclinationAngleSlider.ValueChangedFcn = createCallbackFcn(app, @InclinationAngleSliderValueChanged, true);
            app.InclinationAngleSlider.ValueChangingFcn = createCallbackFcn(app, @InclinationAngleSliderValueChanging, true);
            app.InclinationAngleSlider.BusyAction = 'cancel';
            app.InclinationAngleSlider.Position = [394 130 150 3];

            % Create PhotonshotnoiseSwitchLabel
            app.PhotonshotnoiseSwitchLabel = uilabel(app.FluorophoreTab);
            app.PhotonshotnoiseSwitchLabel.Position = [12 227 102 22];
            app.PhotonshotnoiseSwitchLabel.Text = 'Photon shot noise';

            % Create PhotonshotnoiseSwitch
            app.PhotonshotnoiseSwitch = uiswitch(app.FluorophoreTab, 'slider');
            app.PhotonshotnoiseSwitch.ValueChangedFcn = createCallbackFcn(app, @PhotonshotnoiseSwitchValueChanged, true);
            app.PhotonshotnoiseSwitch.BusyAction = 'cancel';
            app.PhotonshotnoiseSwitch.Position = [160 228 45 20];

            % Create PhotonnumberSpinnerLabel
            app.PhotonnumberSpinnerLabel = uilabel(app.FluorophoreTab);
            app.PhotonnumberSpinnerLabel.Position = [12 290 85 22];
            app.PhotonnumberSpinnerLabel.Text = 'Photon number';

            % Create PhotonnumberSpinner
            app.PhotonnumberSpinner = uispinner(app.FluorophoreTab);
            app.PhotonnumberSpinner.Step = 10000;
            app.PhotonnumberSpinner.ValueChangingFcn = createCallbackFcn(app, @PhotonnumberSpinnerValueChanging, true);
            app.PhotonnumberSpinner.Limits = [0 Inf];
            app.PhotonnumberSpinner.ValueDisplayFormat = '%.0f';
            app.PhotonnumberSpinner.ValueChangedFcn = createCallbackFcn(app, @PhotonnumberSpinnerValueChanged, true);
            app.PhotonnumberSpinner.BusyAction = 'cancel';
            app.PhotonnumberSpinner.Position = [133 290 100 22];
            app.PhotonnumberSpinner.Value = 50000;

            % Create ReducedexcitationSwitchLabel
            app.ReducedexcitationSwitchLabel = uilabel(app.FluorophoreTab);
            app.ReducedexcitationSwitchLabel.Tooltip = {'Consider fluorophore dipole orientation when calculating emission intensity'};
            app.ReducedexcitationSwitchLabel.Position = [12 257 108 22];
            app.ReducedexcitationSwitchLabel.Text = 'Reduced excitation';

            % Create ReducedexcitationSwitch
            app.ReducedexcitationSwitch = uiswitch(app.FluorophoreTab, 'slider');
            app.ReducedexcitationSwitch.ValueChangedFcn = createCallbackFcn(app, @ReducedexcitationSwitchValueChanged, true);
            app.ReducedexcitationSwitch.BusyAction = 'cancel';
            app.ReducedexcitationSwitch.Position = [160 258 45 20];

            % Create ConfiguremultiplefluorophoresButton
            app.ConfiguremultiplefluorophoresButton = uibutton(app.FluorophoreTab, 'push');
            app.ConfiguremultiplefluorophoresButton.ButtonPushedFcn = createCallbackFcn(app, @ConfiguremultiplefluorophoresButtonPushed, true);
            app.ConfiguremultiplefluorophoresButton.BusyAction = 'cancel';
            app.ConfiguremultiplefluorophoresButton.Position = [17 25 182 22];
            app.ConfiguremultiplefluorophoresButton.Text = 'Configure multiple fluorophores';

            % Create zpositionSpinnerLabel
            app.zpositionSpinnerLabel = uilabel(app.FluorophoreTab);
            app.zpositionSpinnerLabel.Position = [12 86 59 22];
            app.zpositionSpinnerLabel.Text = 'z-position';

            % Create zpositionSpinner
            app.zpositionSpinner = uispinner(app.FluorophoreTab);
            app.zpositionSpinner.Step = 10;
            app.zpositionSpinner.ValueChangingFcn = createCallbackFcn(app, @zpositionSpinnerValueChanging, true);
            app.zpositionSpinner.Limits = [0 Inf];
            app.zpositionSpinner.ValueDisplayFormat = '%11.4g nm';
            app.zpositionSpinner.ValueChangedFcn = createCallbackFcn(app, @zpositionSpinnerValueChanged, true);
            app.zpositionSpinner.BusyAction = 'cancel';
            app.zpositionSpinner.Position = [132 86 102 22];

            % Create theta
            app.theta = uieditfield(app.FluorophoreTab, 'numeric');
            app.theta.Limits = [0 90];
            app.theta.RoundFractionalValues = 'on';
            app.theta.ValueDisplayFormat = '%d°';
            app.theta.ValueChangedFcn = createCallbackFcn(app, @thetaValueChanged, true);
            app.theta.Position = [555 118 47 22];

            % Create phi
            app.phi = uieditfield(app.FluorophoreTab, 'numeric');
            app.phi.Limits = [0 360];
            app.phi.RoundFractionalValues = 'on';
            app.phi.ValueDisplayFormat = '%d°';
            app.phi.ValueChangedFcn = createCallbackFcn(app, @phiValueChanged, true);
            app.phi.Position = [556 74 47 22];

            % Create DipolerotationButtonGroup
            app.DipolerotationButtonGroup = uibuttongroup(app.FluorophoreTab);
            app.DipolerotationButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @DipolerotationButtonGroupSelectionChanged, true);
            app.DipolerotationButtonGroup.BorderType = 'none';
            app.DipolerotationButtonGroup.BusyAction = 'cancel';
            app.DipolerotationButtonGroup.Position = [387 153 162 30];

            % Create freelyrotatingButton
            app.freelyrotatingButton = uiradiobutton(app.DipolerotationButtonGroup);
            app.freelyrotatingButton.Tooltip = {'Freely rotating dipole results in isotropic PSF'};
            app.freelyrotatingButton.Text = 'freely rotating';
            app.freelyrotatingButton.Position = [63 5 95 22];

            % Create fixedButton
            app.fixedButton = uiradiobutton(app.DipolerotationButtonGroup);
            app.fixedButton.Text = 'fixed';
            app.fixedButton.Position = [7 5 47 22];
            app.fixedButton.Value = true;

            % Create DipoleorientationLabel
            app.DipoleorientationLabel = uilabel(app.FluorophoreTab);
            app.DipoleorientationLabel.Position = [273 157 98 22];
            app.DipoleorientationLabel.Text = 'Dipole orientation';

            % Create EmissionwavelengthSpinnerLabel
            app.EmissionwavelengthSpinnerLabel = uilabel(app.FluorophoreTab);
            app.EmissionwavelengthSpinnerLabel.Position = [12 321 118 22];
            app.EmissionwavelengthSpinnerLabel.Text = 'Emission wavelength';

            % Create EmissionWavelengthSpinner
            app.EmissionWavelengthSpinner = uispinner(app.FluorophoreTab);
            app.EmissionWavelengthSpinner.Step = 10;
            app.EmissionWavelengthSpinner.Limits = [200 Inf];
            app.EmissionWavelengthSpinner.ValueDisplayFormat = '%11.4g nm';
            app.EmissionWavelengthSpinner.ValueChangedFcn = createCallbackFcn(app, @EmissionWavelengthSpinnerValueChanged, true);
            app.EmissionWavelengthSpinner.BusyAction = 'cancel';
            app.EmissionWavelengthSpinner.Position = [133 321 100 22];
            app.EmissionWavelengthSpinner.Value = 680;

            % Create ypositionSpinnerLabel
            app.ypositionSpinnerLabel = uilabel(app.FluorophoreTab);
            app.ypositionSpinnerLabel.Position = [12 121 59 22];
            app.ypositionSpinnerLabel.Text = 'y-position';

            % Create ypositionSpinner
            app.ypositionSpinner = uispinner(app.FluorophoreTab);
            app.ypositionSpinner.Step = 10;
            app.ypositionSpinner.ValueChangingFcn = createCallbackFcn(app, @ypositionSpinnerValueChanging, true);
            app.ypositionSpinner.ValueDisplayFormat = '%11.4g nm';
            app.ypositionSpinner.ValueChangedFcn = createCallbackFcn(app, @ypositionSpinnerValueChanged, true);
            app.ypositionSpinner.BusyAction = 'cancel';
            app.ypositionSpinner.Position = [132 121 102 22];

            % Create xpositionSpinnerLabel
            app.xpositionSpinnerLabel = uilabel(app.FluorophoreTab);
            app.xpositionSpinnerLabel.Position = [12 157 59 22];
            app.xpositionSpinnerLabel.Text = 'x-position';

            % Create xpositionSpinner
            app.xpositionSpinner = uispinner(app.FluorophoreTab);
            app.xpositionSpinner.Step = 10;
            app.xpositionSpinner.ValueChangingFcn = createCallbackFcn(app, @xpositionSpinnerValueChanging, true);
            app.xpositionSpinner.ValueDisplayFormat = '%11.4g nm';
            app.xpositionSpinner.ValueChangedFcn = createCallbackFcn(app, @xpositionSpinnerValueChanged, true);
            app.xpositionSpinner.BusyAction = 'cancel';
            app.xpositionSpinner.Position = [132 157 102 22];

            % Create PositionLabel
            app.PositionLabel = uilabel(app.FluorophoreTab);
            app.PositionLabel.FontSize = 13;
            app.PositionLabel.FontWeight = 'bold';
            app.PositionLabel.Position = [13 187 82 22];
            app.PositionLabel.Text = 'Position';

            % Create OrientationLabel
            app.OrientationLabel = uilabel(app.FluorophoreTab);
            app.OrientationLabel.FontSize = 13;
            app.OrientationLabel.FontWeight = 'bold';
            app.OrientationLabel.Position = [273 187 82 22];
            app.OrientationLabel.Text = 'Orientation';

            % Create MicroscopeRITab
            app.MicroscopeRITab = uitab(app.TabGroup);
            app.MicroscopeRITab.Title = 'Microscope & RI';

            % Create ObjectiveandtubelensLabel
            app.ObjectiveandtubelensLabel = uilabel(app.MicroscopeRITab);
            app.ObjectiveandtubelensLabel.FontSize = 13;
            app.ObjectiveandtubelensLabel.FontWeight = 'bold';
            app.ObjectiveandtubelensLabel.Position = [16 345 151 22];
            app.ObjectiveandtubelensLabel.Text = 'Objective and tube lens';

            % Create ObjectiveNaLabel
            app.ObjectiveNaLabel = uilabel(app.MicroscopeRITab);
            app.ObjectiveNaLabel.Position = [16 317 85 22];
            app.ObjectiveNaLabel.Text = 'Objective NA';

            % Create DefocusSliderLabel
            app.DefocusSliderLabel = uilabel(app.MicroscopeRITab);
            app.DefocusSliderLabel.Position = [16 246 50 22];
            app.DefocusSliderLabel.Text = 'Defocus';

            % Create DefocusSlider
            app.DefocusSlider = uislider(app.MicroscopeRITab);
            app.DefocusSlider.Limits = [-1000 1000];
            app.DefocusSlider.MajorTicks = [-1000 0 1000];
            app.DefocusSlider.ValueChangedFcn = createCallbackFcn(app, @DefocusSliderValueChanged, true);
            app.DefocusSlider.ValueChangingFcn = createCallbackFcn(app, @DefocusSliderValueChanging, true);
            app.DefocusSlider.BusyAction = 'cancel';
            app.DefocusSlider.Position = [129 255 178 3];

            % Create ObjectiveNALabel
            app.ObjectiveNALabel = uilabel(app.MicroscopeRITab);
            app.ObjectiveNALabel.Position = [16 317 85 22];
            app.ObjectiveNALabel.Text = '';

            % Create ObjectiveNaSpinner
            app.ObjectiveNaSpinner = uispinner(app.MicroscopeRITab);
            app.ObjectiveNaSpinner.Step = 0.1;
            app.ObjectiveNaSpinner.ValueChangingFcn = createCallbackFcn(app, @ObjectiveNaSpinnerValueChanging, true);
            app.ObjectiveNaSpinner.Limits = [0.1 2];
            app.ObjectiveNaSpinner.ValueChangedFcn = createCallbackFcn(app, @ObjectiveNaSpinnerValueChanged, true);
            app.ObjectiveNaSpinner.BusyAction = 'cancel';
            app.ObjectiveNaSpinner.Position = [152 317 77 22];
            app.ObjectiveNaSpinner.Value = 0.7;

            % Create FocusButton
            app.FocusButton = uibutton(app.MicroscopeRITab, 'push');
            app.FocusButton.ButtonPushedFcn = createCallbackFcn(app, @FocusButtonPushed, true);
            app.FocusButton.BusyAction = 'cancel';
            app.FocusButton.Tooltip = {'Set the defocus value to 0nm'};
            app.FocusButton.Position = [422 245 71 23];
            app.FocusButton.Text = 'Focus';

            % Create defocus
            app.defocus = uieditfield(app.MicroscopeRITab, 'numeric');
            app.defocus.Limits = [-1000 1000];
            app.defocus.RoundFractionalValues = 'on';
            app.defocus.ValueDisplayFormat = '%11.4g nm';
            app.defocus.ValueChangedFcn = createCallbackFcn(app, @defocusValueChanged, true);
            app.defocus.Position = [330 246 69 22];

            % Create RefractiveindexLabel
            app.RefractiveindexLabel = uilabel(app.MicroscopeRITab);
            app.RefractiveindexLabel.FontSize = 13;
            app.RefractiveindexLabel.FontWeight = 'bold';
            app.RefractiveindexLabel.Position = [12 178 116 22];
            app.RefractiveindexLabel.Text = 'Refractive indices';

            % Create RefractiveIndexSpecimenLabel
            app.RefractiveIndexSpecimenLabel = uilabel(app.MicroscopeRITab);
            app.RefractiveIndexSpecimenLabel.Position = [12 123 85 22];
            app.RefractiveIndexSpecimenLabel.Text = 'RI Specimen';

            % Create RefractiveIndexSpecimenSpinner
            app.RefractiveIndexSpecimenSpinner = uispinner(app.MicroscopeRITab);
            app.RefractiveIndexSpecimenSpinner.Step = 0.1;
            app.RefractiveIndexSpecimenSpinner.Limits = [1 5];
            app.RefractiveIndexSpecimenSpinner.ValueChangedFcn = createCallbackFcn(app, @RefractiveIndexSpecimenSpinnerValueChanged, true);
            app.RefractiveIndexSpecimenSpinner.BusyAction = 'cancel';
            app.RefractiveIndexSpecimenSpinner.Position = [169 123 77 22];
            app.RefractiveIndexSpecimenSpinner.Value = 1.33;

            % Create RefractiveIndexIntermediateLayerLabel
            app.RefractiveIndexIntermediateLayerLabel = uilabel(app.MicroscopeRITab);
            app.RefractiveIndexIntermediateLayerLabel.Visible = 'off';
            app.RefractiveIndexIntermediateLayerLabel.Position = [12 67 116 22];
            app.RefractiveIndexIntermediateLayerLabel.Text = 'RI Intermediate layer';

            % Create RefractiveIndexIntermediateLayerSpinner
            app.RefractiveIndexIntermediateLayerSpinner = uispinner(app.MicroscopeRITab);
            app.RefractiveIndexIntermediateLayerSpinner.Step = 0.1;
            app.RefractiveIndexIntermediateLayerSpinner.ValueChangingFcn = createCallbackFcn(app, @RefractiveIndexIntermediateLayerSpinnerValueChanging, true);
            app.RefractiveIndexIntermediateLayerSpinner.Limits = [1 5];
            app.RefractiveIndexIntermediateLayerSpinner.ValueChangedFcn = createCallbackFcn(app, @RefractiveIndexIntermediateLayerSpinnerValueChanged, true);
            app.RefractiveIndexIntermediateLayerSpinner.BusyAction = 'cancel';
            app.RefractiveIndexIntermediateLayerSpinner.Visible = 'off';
            app.RefractiveIndexIntermediateLayerSpinner.Position = [169 67 77 22];
            app.RefractiveIndexIntermediateLayerSpinner.Value = 1.46;

            % Create RefractiveIndexImmersionMediumLabel
            app.RefractiveIndexImmersionMediumLabel = uilabel(app.MicroscopeRITab);
            app.RefractiveIndexImmersionMediumLabel.Position = [12 95 122 22];
            app.RefractiveIndexImmersionMediumLabel.Text = 'RI Immersion medium';

            % Create RefractiveIndexImmersionMediumSpinner
            app.RefractiveIndexImmersionMediumSpinner = uispinner(app.MicroscopeRITab);
            app.RefractiveIndexImmersionMediumSpinner.Step = 0.1;
            app.RefractiveIndexImmersionMediumSpinner.ValueChangingFcn = createCallbackFcn(app, @RefractiveIndexImmersionMediumSpinnerValueChanging, true);
            app.RefractiveIndexImmersionMediumSpinner.Limits = [1 5];
            app.RefractiveIndexImmersionMediumSpinner.ValueChangedFcn = createCallbackFcn(app, @RefractiveIndexImmersionMediumSpinnerValueChanged, true);
            app.RefractiveIndexImmersionMediumSpinner.BusyAction = 'cancel';
            app.RefractiveIndexImmersionMediumSpinner.Position = [169 95 77 22];
            app.RefractiveIndexImmersionMediumSpinner.Value = 1;

            % Create MagnificationLabel
            app.MagnificationLabel = uilabel(app.MicroscopeRITab);
            app.MagnificationLabel.Tooltip = {'Magnification = Focal length tube lens / Focal length objective'};
            app.MagnificationLabel.Position = [16 287 85 22];
            app.MagnificationLabel.Text = 'Magnification';

            % Create MagnificationSpinner
            app.MagnificationSpinner = uispinner(app.MicroscopeRITab);
            app.MagnificationSpinner.Step = 10;
            app.MagnificationSpinner.Limits = [1 Inf];
            app.MagnificationSpinner.ValueDisplayFormat = '%11.4gx';
            app.MagnificationSpinner.ValueChangedFcn = createCallbackFcn(app, @MagnificationSpinnerValueChanged, true);
            app.MagnificationSpinner.BusyAction = 'cancel';
            app.MagnificationSpinner.Position = [152 287 77 22];
            app.MagnificationSpinner.Value = 60;

            % Create AddintermediatelayerCheckBox
            app.AddintermediatelayerCheckBox = uicheckbox(app.MicroscopeRITab);
            app.AddintermediatelayerCheckBox.ValueChangedFcn = createCallbackFcn(app, @AddintermediatelayerCheckBoxValueChanged, true);
            app.AddintermediatelayerCheckBox.Text = 'Add intermediate layer';
            app.AddintermediatelayerCheckBox.Position = [12 150 142 22];

            % Create IntermediateLayerThicknessLabel
            app.IntermediateLayerThicknessLabel = uilabel(app.MicroscopeRITab);
            app.IntermediateLayerThicknessLabel.Visible = 'off';
            app.IntermediateLayerThicknessLabel.Position = [12 39 154 22];
            app.IntermediateLayerThicknessLabel.Text = 'Intermediate layer thickness';

            % Create IntermediateLayerThicknessSpinner
            app.IntermediateLayerThicknessSpinner = uispinner(app.MicroscopeRITab);
            app.IntermediateLayerThicknessSpinner.Step = 0.1;
            app.IntermediateLayerThicknessSpinner.Limits = [0 Inf];
            app.IntermediateLayerThicknessSpinner.ValueDisplayFormat = '%11.4g µm';
            app.IntermediateLayerThicknessSpinner.ValueChangedFcn = createCallbackFcn(app, @IntermediateLayerThicknessSpinnerValueChanged, true);
            app.IntermediateLayerThicknessSpinner.BusyAction = 'cancel';
            app.IntermediateLayerThicknessSpinner.Visible = 'off';
            app.IntermediateLayerThicknessSpinner.Position = [169 39 77 22];

            % Create CameraLabel
            app.CameraLabel = uilabel(app.MicroscopeRITab);
            app.CameraLabel.FontSize = 13;
            app.CameraLabel.FontWeight = 'bold';
            app.CameraLabel.Position = [267 174 53 22];
            app.CameraLabel.Text = 'Camera';

            % Create BackgroundnoisestdSpinnerLabel
            app.BackgroundnoisestdSpinnerLabel = uilabel(app.MicroscopeRITab);
            app.BackgroundnoisestdSpinnerLabel.Position = [267 90 129 22];
            app.BackgroundnoisestdSpinnerLabel.Text = 'Background noise (std)';

            % Create BackgroundnoisestdSpinner
            app.BackgroundnoisestdSpinner = uispinner(app.MicroscopeRITab);
            app.BackgroundnoisestdSpinner.Step = 10;
            app.BackgroundnoisestdSpinner.ValueChangingFcn = createCallbackFcn(app, @BackgroundnoisestdSpinnerValueChanging, true);
            app.BackgroundnoisestdSpinner.Limits = [0 Inf];
            app.BackgroundnoisestdSpinner.ValueDisplayFormat = '%.0f';
            app.BackgroundnoisestdSpinner.ValueChangedFcn = createCallbackFcn(app, @BackgroundnoisestdSpinnerValueChanged, true);
            app.BackgroundnoisestdSpinner.BusyAction = 'cancel';
            app.BackgroundnoisestdSpinner.Position = [403 90 90 22];

            % Create PixelsizephysicalSpinnerLabel
            app.PixelsizephysicalSpinnerLabel = uilabel(app.MicroscopeRITab);
            app.PixelsizephysicalSpinnerLabel.Position = [267 145 110 22];
            app.PixelsizephysicalSpinnerLabel.Text = 'Pixel size (physical)';

            % Create PixelsizePhysicalSpinner
            app.PixelsizePhysicalSpinner = uispinner(app.MicroscopeRITab);
            app.PixelsizePhysicalSpinner.Step = 0.1;
            app.PixelsizePhysicalSpinner.ValueChangingFcn = createCallbackFcn(app, @PixelsizePhysicalSpinnerValueChanging, true);
            app.PixelsizePhysicalSpinner.Limits = [0.1 Inf];
            app.PixelsizePhysicalSpinner.ValueDisplayFormat = '%11.4g µm';
            app.PixelsizePhysicalSpinner.ValueChangedFcn = createCallbackFcn(app, @PixelsizePhysicalSpinnerValueChanged, true);
            app.PixelsizePhysicalSpinner.BusyAction = 'cancel';
            app.PixelsizePhysicalSpinner.Position = [403 145 90 22];
            app.PixelsizePhysicalSpinner.Value = 6;

            % Create PixelsizeobjectspaceLabel
            app.PixelsizeobjectspaceLabel = uilabel(app.MicroscopeRITab);
            app.PixelsizeobjectspaceLabel.Position = [267 116 134 22];
            app.PixelsizeobjectspaceLabel.Text = 'Pixel size (object space)';

            % Create PixelsizeObjectSpaceSpinner
            app.PixelsizeObjectSpaceSpinner = uispinner(app.MicroscopeRITab);
            app.PixelsizeObjectSpaceSpinner.Step = 10;
            app.PixelsizeObjectSpaceSpinner.Limits = [10 Inf];
            app.PixelsizeObjectSpaceSpinner.ValueDisplayFormat = '%11.4g nm';
            app.PixelsizeObjectSpaceSpinner.ValueChangedFcn = createCallbackFcn(app, @PixelsizeObjectSpaceSpinnerValueChanged, true);
            app.PixelsizeObjectSpaceSpinner.BusyAction = 'cancel';
            app.PixelsizeObjectSpaceSpinner.Position = [403 116 90 22];
            app.PixelsizeObjectSpaceSpinner.Value = 100;

            % Create TubelensfocallengthLabel
            app.TubelensfocallengthLabel = uilabel(app.MicroscopeRITab);
            app.TubelensfocallengthLabel.Position = [269 288 122 22];
            app.TubelensfocallengthLabel.Text = 'Tube lens focal length';

            % Create TubeLensFocalLengthSpinner
            app.TubeLensFocalLengthSpinner = uispinner(app.MicroscopeRITab);
            app.TubeLensFocalLengthSpinner.Limits = [0.1 Inf];
            app.TubeLensFocalLengthSpinner.ValueDisplayFormat = '%11.4g mm';
            app.TubeLensFocalLengthSpinner.ValueChangedFcn = createCallbackFcn(app, @TubeLensFocalLengthSpinnerValueChanged, true);
            app.TubeLensFocalLengthSpinner.BusyAction = 'cancel';
            app.TubeLensFocalLengthSpinner.Position = [405 288 90 22];
            app.TubeLensFocalLengthSpinner.Value = 180;

            % Create ObjectivefocallengthLabel
            app.ObjectivefocallengthLabel = uilabel(app.MicroscopeRITab);
            app.ObjectivefocallengthLabel.Position = [269 317 120 22];
            app.ObjectivefocallengthLabel.Text = 'Objective focal length';

            % Create ObjectiveFocalLengthSpinner
            app.ObjectiveFocalLengthSpinner = uispinner(app.MicroscopeRITab);
            app.ObjectiveFocalLengthSpinner.Limits = [0.1 Inf];
            app.ObjectiveFocalLengthSpinner.ValueDisplayFormat = '%11.4g mm';
            app.ObjectiveFocalLengthSpinner.ValueChangedFcn = createCallbackFcn(app, @ObjectiveFocalLengthSpinnerValueChanged, true);
            app.ObjectiveFocalLengthSpinner.BusyAction = 'cancel';
            app.ObjectiveFocalLengthSpinner.Position = [405 317 90 22];
            app.ObjectiveFocalLengthSpinner.Value = 3;

            % Create AberrationsTab
            app.AberrationsTab = uitab(app.TabGroup);
            app.AberrationsTab.Title = 'Aberrations';

            % Create SpecifyZernikeaberrationsButtonGroup
            app.SpecifyZernikeaberrationsButtonGroup = uibuttongroup(app.AberrationsTab);
            app.SpecifyZernikeaberrationsButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @SpecifyZernikeaberrationsButtonGroupSelectionChanged, true);
            app.SpecifyZernikeaberrationsButtonGroup.Title = 'Specify Zernike aberrations';
            app.SpecifyZernikeaberrationsButtonGroup.Position = [12 112 588 233];

            % Create ZernikeCommonAberrationsButton
            app.ZernikeCommonAberrationsButton = uiradiobutton(app.SpecifyZernikeaberrationsButtonGroup);
            app.ZernikeCommonAberrationsButton.Text = 'Select common aberrations';
            app.ZernikeCommonAberrationsButton.Position = [11 181 168 22];
            app.ZernikeCommonAberrationsButton.Value = true;

            % Create ZernikeInputVectorButton
            app.ZernikeInputVectorButton = uiradiobutton(app.SpecifyZernikeaberrationsButtonGroup);
            app.ZernikeInputVectorButton.Tooltip = {''};
            app.ZernikeInputVectorButton.Text = 'Specify Zernike indices (Noll) and coefficients via input vector';
            app.ZernikeInputVectorButton.Position = [11 44 356 22];

            % Create LoadZernikeIndicesWeightsButton
            app.LoadZernikeIndicesWeightsButton = uibutton(app.SpecifyZernikeaberrationsButtonGroup, 'push');
            app.LoadZernikeIndicesWeightsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadZernikeIndicesWeightsButtonPushed, true);
            app.LoadZernikeIndicesWeightsButton.Enable = 'off';
            app.LoadZernikeIndicesWeightsButton.Tooltip = {'Load Zernike indices and weights from file'};
            app.LoadZernikeIndicesWeightsButton.Position = [455 16 53 23];
            app.LoadZernikeIndicesWeightsButton.Text = 'Load';

            % Create IndicesWeightsEditField
            app.IndicesWeightsEditField = uieditfield(app.SpecifyZernikeaberrationsButtonGroup, 'text');
            app.IndicesWeightsEditField.ValueChangedFcn = createCallbackFcn(app, @IndicesWeightsEditFieldValueChanged, true);
            app.IndicesWeightsEditField.Enable = 'off';
            app.IndicesWeightsEditField.Tooltip = {'Specify Zernike Noll indices and coefficients in the format index:weight separated by semicolons. Zernike coefficients should be specified in units of the wavelength.'};
            app.IndicesWeightsEditField.Placeholder = '1:0.0;2:0.0;3:0.0;4:0.0';
            app.IndicesWeightsEditField.Position = [30 16 403 22];

            % Create SaveZernikeIndicesWeightsButton
            app.SaveZernikeIndicesWeightsButton = uibutton(app.SpecifyZernikeaberrationsButtonGroup, 'push');
            app.SaveZernikeIndicesWeightsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveZernikeIndicesWeightsButtonPushed, true);
            app.SaveZernikeIndicesWeightsButton.Enable = 'off';
            app.SaveZernikeIndicesWeightsButton.Tooltip = {'Load Zernike indices and weights from file'};
            app.SaveZernikeIndicesWeightsButton.Position = [514 16 53 23];
            app.SaveZernikeIndicesWeightsButton.Text = 'Save';

            % Create ZernikeaberrationsLabel
            app.ZernikeaberrationsLabel = uilabel(app.AberrationsTab);
            app.ZernikeaberrationsLabel.FontSize = 13;
            app.ZernikeaberrationsLabel.FontWeight = 'bold';
            app.ZernikeaberrationsLabel.Position = [16 352 126 22];
            app.ZernikeaberrationsLabel.Text = 'Zernike aberrations';

            % Create ZernikeAberrationsShowplotCheckBox
            app.ZernikeAberrationsShowplotCheckBox = uicheckbox(app.AberrationsTab);
            app.ZernikeAberrationsShowplotCheckBox.ValueChangedFcn = createCallbackFcn(app, @ZernikeAberrationsShowplotCheckBoxValueChanged, true);
            app.ZernikeAberrationsShowplotCheckBox.Text = 'Show plot';
            app.ZernikeAberrationsShowplotCheckBox.Position = [166 352 76 22];

            % Create FitZernikeaberrationsLabel
            app.FitZernikeaberrationsLabel = uilabel(app.AberrationsTab);
            app.FitZernikeaberrationsLabel.FontSize = 13;
            app.FitZernikeaberrationsLabel.FontWeight = 'bold';
            app.FitZernikeaberrationsLabel.Position = [16 76 145 22];
            app.FitZernikeaberrationsLabel.Text = 'Fit Zernike aberrations';

            % Create FitZernikeLabel
            app.FitZernikeLabel = uilabel(app.AberrationsTab);
            app.FitZernikeLabel.Position = [16 46 291 22];
            app.FitZernikeLabel.Text = 'Estimate Zernike coefficients from experimental data:';

            % Create FitZernikeButton
            app.FitZernikeButton = uibutton(app.AberrationsTab, 'push');
            app.FitZernikeButton.ButtonPushedFcn = createCallbackFcn(app, @FitZernikeButtonPushed, true);
            app.FitZernikeButton.Position = [311 45 100 23];
            app.FitZernikeButton.Text = 'Fit Zernike';

            % Create HorizontalComaCheckBox
            app.HorizontalComaCheckBox = uicheckbox(app.AberrationsTab);
            app.HorizontalComaCheckBox.ValueChangedFcn = createCallbackFcn(app, @HorizontalComaCheckBoxValueChanged, true);
            app.HorizontalComaCheckBox.Text = 'Horizontal coma';
            app.HorizontalComaCheckBox.Position = [337 197 128 22];

            % Create VerticalComaCheckBox
            app.VerticalComaCheckBox = uicheckbox(app.AberrationsTab);
            app.VerticalComaCheckBox.ValueChangedFcn = createCallbackFcn(app, @VerticalComaCheckBoxValueChanged, true);
            app.VerticalComaCheckBox.Text = 'Vertical coma';
            app.VerticalComaCheckBox.Position = [337 220 128 22];

            % Create VerticalAstigmatismCheckBox
            app.VerticalAstigmatismCheckBox = uicheckbox(app.AberrationsTab);
            app.VerticalAstigmatismCheckBox.ValueChangedFcn = createCallbackFcn(app, @VerticalAstigmatismCheckBoxValueChanged, true);
            app.VerticalAstigmatismCheckBox.Text = 'Vertical astigmatism';
            app.VerticalAstigmatismCheckBox.Position = [337 243 128 22];

            % Create ObliqueAstigmatismCheckBox
            app.ObliqueAstigmatismCheckBox = uicheckbox(app.AberrationsTab);
            app.ObliqueAstigmatismCheckBox.ValueChangedFcn = createCallbackFcn(app, @ObliqueAstigmatismCheckBoxValueChanged, true);
            app.ObliqueAstigmatismCheckBox.Text = 'Oblique astigmatism';
            app.ObliqueAstigmatismCheckBox.Position = [337 266 128 22];

            % Create PrimarySphericalCheckBox
            app.PrimarySphericalCheckBox = uicheckbox(app.AberrationsTab);
            app.PrimarySphericalCheckBox.ValueChangedFcn = createCallbackFcn(app, @PrimarySphericalCheckBoxValueChanged, true);
            app.PrimarySphericalCheckBox.Text = 'Primary spherical';
            app.PrimarySphericalCheckBox.Position = [42 195 128 22];

            % Create DefocusCheckBox
            app.DefocusCheckBox = uicheckbox(app.AberrationsTab);
            app.DefocusCheckBox.ValueChangedFcn = createCallbackFcn(app, @DefocusCheckBoxValueChanged, true);
            app.DefocusCheckBox.Text = 'Defocus';
            app.DefocusCheckBox.Position = [42 219 128 22];

            % Create HorizontalTiltCheckBox
            app.HorizontalTiltCheckBox = uicheckbox(app.AberrationsTab);
            app.HorizontalTiltCheckBox.ValueChangedFcn = createCallbackFcn(app, @HorizontalTiltCheckBoxValueChanged, true);
            app.HorizontalTiltCheckBox.Text = 'Horizontal tilt';
            app.HorizontalTiltCheckBox.Position = [42 243 128 22];

            % Create WeightVerticalTiltSpinnerLabel
            app.WeightVerticalTiltSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightVerticalTiltSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightVerticalTiltSpinnerLabel.Enable = 'off';
            app.WeightVerticalTiltSpinnerLabel.Position = [171 267 25 22];
            app.WeightVerticalTiltSpinnerLabel.Text = 'w =';

            % Create WeightVerticalTiltSpinner
            app.WeightVerticalTiltSpinner = uispinner(app.AberrationsTab);
            app.WeightVerticalTiltSpinner.Step = 0.01;
            app.WeightVerticalTiltSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightVerticalTiltSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightVerticalTiltSpinnerValueChanged, true);
            app.WeightVerticalTiltSpinner.BusyAction = 'cancel';
            app.WeightVerticalTiltSpinner.Enable = 'off';
            app.WeightVerticalTiltSpinner.Position = [201 267 67 22];

            % Create VerticalTiltCheckBox
            app.VerticalTiltCheckBox = uicheckbox(app.AberrationsTab);
            app.VerticalTiltCheckBox.ValueChangedFcn = createCallbackFcn(app, @VerticalTiltCheckBoxValueChanged, true);
            app.VerticalTiltCheckBox.Text = 'Vertical tilt';
            app.VerticalTiltCheckBox.Position = [42 267 128 22];

            % Create WeightHorizontalTiltSpinnerLabel
            app.WeightHorizontalTiltSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightHorizontalTiltSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightHorizontalTiltSpinnerLabel.Enable = 'off';
            app.WeightHorizontalTiltSpinnerLabel.Position = [171 243 25 22];
            app.WeightHorizontalTiltSpinnerLabel.Text = 'w =';

            % Create WeightHorizontalTiltSpinner
            app.WeightHorizontalTiltSpinner = uispinner(app.AberrationsTab);
            app.WeightHorizontalTiltSpinner.Step = 0.01;
            app.WeightHorizontalTiltSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightHorizontalTiltSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightHorizontalTiltSpinnerValueChanged, true);
            app.WeightHorizontalTiltSpinner.BusyAction = 'cancel';
            app.WeightHorizontalTiltSpinner.Enable = 'off';
            app.WeightHorizontalTiltSpinner.Position = [201 243 67 22];

            % Create WeightDefocusSpinnerLabel
            app.WeightDefocusSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightDefocusSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightDefocusSpinnerLabel.Enable = 'off';
            app.WeightDefocusSpinnerLabel.Position = [171 219 25 22];
            app.WeightDefocusSpinnerLabel.Text = 'w =';

            % Create WeightDefocusSpinner
            app.WeightDefocusSpinner = uispinner(app.AberrationsTab);
            app.WeightDefocusSpinner.Step = 0.01;
            app.WeightDefocusSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightDefocusSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightDefocusSpinnerValueChanged, true);
            app.WeightDefocusSpinner.BusyAction = 'cancel';
            app.WeightDefocusSpinner.Enable = 'off';
            app.WeightDefocusSpinner.Position = [201 219 67 22];

            % Create WeightPrimarySphericalSpinnerLabel
            app.WeightPrimarySphericalSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightPrimarySphericalSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightPrimarySphericalSpinnerLabel.Enable = 'off';
            app.WeightPrimarySphericalSpinnerLabel.Position = [171 195 25 22];
            app.WeightPrimarySphericalSpinnerLabel.Text = 'w =';

            % Create WeightPrimarySphericalSpinner
            app.WeightPrimarySphericalSpinner = uispinner(app.AberrationsTab);
            app.WeightPrimarySphericalSpinner.Step = 0.01;
            app.WeightPrimarySphericalSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightPrimarySphericalSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightPrimarySphericalSpinnerValueChanged, true);
            app.WeightPrimarySphericalSpinner.BusyAction = 'cancel';
            app.WeightPrimarySphericalSpinner.Enable = 'off';
            app.WeightPrimarySphericalSpinner.Position = [201 195 67 22];

            % Create WeightObliqueAstigmatismSpinnerLabel
            app.WeightObliqueAstigmatismSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightObliqueAstigmatismSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightObliqueAstigmatismSpinnerLabel.Enable = 'off';
            app.WeightObliqueAstigmatismSpinnerLabel.Position = [485 266 25 22];
            app.WeightObliqueAstigmatismSpinnerLabel.Text = 'w =';

            % Create WeightObliqueAstigmatismSpinner
            app.WeightObliqueAstigmatismSpinner = uispinner(app.AberrationsTab);
            app.WeightObliqueAstigmatismSpinner.Step = 0.01;
            app.WeightObliqueAstigmatismSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightObliqueAstigmatismSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightObliqueAstigmatismSpinnerValueChanged, true);
            app.WeightObliqueAstigmatismSpinner.BusyAction = 'cancel';
            app.WeightObliqueAstigmatismSpinner.Enable = 'off';
            app.WeightObliqueAstigmatismSpinner.Position = [515 266 67 22];

            % Create WeightVerticalAstigmatismSpinnerLabel
            app.WeightVerticalAstigmatismSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightVerticalAstigmatismSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightVerticalAstigmatismSpinnerLabel.Enable = 'off';
            app.WeightVerticalAstigmatismSpinnerLabel.Position = [485 243 25 22];
            app.WeightVerticalAstigmatismSpinnerLabel.Text = 'w =';

            % Create WeightVerticalAstigmatismSpinner
            app.WeightVerticalAstigmatismSpinner = uispinner(app.AberrationsTab);
            app.WeightVerticalAstigmatismSpinner.Step = 0.01;
            app.WeightVerticalAstigmatismSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightVerticalAstigmatismSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightVerticalAstigmatismSpinnerValueChanged, true);
            app.WeightVerticalAstigmatismSpinner.BusyAction = 'cancel';
            app.WeightVerticalAstigmatismSpinner.Enable = 'off';
            app.WeightVerticalAstigmatismSpinner.Position = [515 243 67 22];

            % Create WeightVerticalComaSpinnerLabel
            app.WeightVerticalComaSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightVerticalComaSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightVerticalComaSpinnerLabel.Enable = 'off';
            app.WeightVerticalComaSpinnerLabel.Position = [485 220 25 22];
            app.WeightVerticalComaSpinnerLabel.Text = 'w =';

            % Create WeightVerticalComaSpinner
            app.WeightVerticalComaSpinner = uispinner(app.AberrationsTab);
            app.WeightVerticalComaSpinner.Step = 0.01;
            app.WeightVerticalComaSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightVerticalComaSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightVerticalComaSpinnerValueChanged, true);
            app.WeightVerticalComaSpinner.BusyAction = 'cancel';
            app.WeightVerticalComaSpinner.Enable = 'off';
            app.WeightVerticalComaSpinner.Position = [515 220 67 22];

            % Create WeightHorizontalComaSpinnerLabel
            app.WeightHorizontalComaSpinnerLabel = uilabel(app.AberrationsTab);
            app.WeightHorizontalComaSpinnerLabel.HorizontalAlignment = 'right';
            app.WeightHorizontalComaSpinnerLabel.Enable = 'off';
            app.WeightHorizontalComaSpinnerLabel.Position = [485 197 25 22];
            app.WeightHorizontalComaSpinnerLabel.Text = 'w =';

            % Create WeightHorizontalComaSpinner
            app.WeightHorizontalComaSpinner = uispinner(app.AberrationsTab);
            app.WeightHorizontalComaSpinner.Step = 0.01;
            app.WeightHorizontalComaSpinner.ValueDisplayFormat = '%11.4g \x03bb';
            app.WeightHorizontalComaSpinner.ValueChangedFcn = createCallbackFcn(app, @WeightHorizontalComaSpinnerValueChanged, true);
            app.WeightHorizontalComaSpinner.BusyAction = 'cancel';
            app.WeightHorizontalComaSpinner.Enable = 'off';
            app.WeightHorizontalComaSpinner.Position = [515 197 67 22];

            % Create PhasemaskTab
            app.PhasemaskTab = uitab(app.TabGroup);
            app.PhasemaskTab.Title = 'Phase mask';

            % Create PhasemaskLabel
            app.PhasemaskLabel = uilabel(app.PhasemaskTab);
            app.PhasemaskLabel.FontSize = 13;
            app.PhasemaskLabel.FontWeight = 'bold';
            app.PhasemaskLabel.Position = [16 352 81 22];
            app.PhasemaskLabel.Text = 'Phase mask';

            % Create BFPmanipulationDropDownLabel
            app.BFPmanipulationDropDownLabel = uilabel(app.PhasemaskTab);
            app.BFPmanipulationDropDownLabel.Position = [16 315 101 22];
            app.BFPmanipulationDropDownLabel.Text = 'BFP manipulation';

            % Create BFPmanipulationDropDown
            app.BFPmanipulationDropDown = uidropdown(app.PhasemaskTab);
            app.BFPmanipulationDropDown.Items = {'None', 'Astigmatism', 'Vortex', 'Sector', 'Opposing Sectors', 'Pyramid', 'Double Helix', 'Custom'};
            app.BFPmanipulationDropDown.ItemsData = {'none', 'Astigmatism', 'Vortex', 'Sector', 'OpposingSectors', 'Pyramid', 'DoubleHelix', 'Custom'};
            app.BFPmanipulationDropDown.ValueChangedFcn = createCallbackFcn(app, @BFPmanipulationDropDownValueChanged, true);
            app.BFPmanipulationDropDown.BusyAction = 'cancel';
            app.BFPmanipulationDropDown.Position = [136 315 150 22];
            app.BFPmanipulationDropDown.Value = 'none';

            % Create InnerringradiusSliderLabel
            app.InnerringradiusSliderLabel = uilabel(app.PhasemaskTab);
            app.InnerringradiusSliderLabel.Visible = 'off';
            app.InnerringradiusSliderLabel.Position = [16 239 92 22];
            app.InnerringradiusSliderLabel.Text = 'Inner ring radius';

            % Create InnerringradiusSlider
            app.InnerringradiusSlider = uislider(app.PhasemaskTab);
            app.InnerringradiusSlider.Limits = [0 1];
            app.InnerringradiusSlider.ValueChangedFcn = createCallbackFcn(app, @InnerringradiusSliderValueChanged, true);
            app.InnerringradiusSlider.ValueChangingFcn = createCallbackFcn(app, @InnerringradiusSliderValueChanging, true);
            app.InnerringradiusSlider.BusyAction = 'cancel';
            app.InnerringradiusSlider.Visible = 'off';
            app.InnerringradiusSlider.Position = [136 248 150 3];

            % Create RotatephasemaskSliderLabel
            app.RotatephasemaskSliderLabel = uilabel(app.PhasemaskTab);
            app.RotatephasemaskSliderLabel.Visible = 'off';
            app.RotatephasemaskSliderLabel.Position = [16 186 109 22];
            app.RotatephasemaskSliderLabel.Text = 'Rotate phase mask';

            % Create RotatephasemaskSlider
            app.RotatephasemaskSlider = uislider(app.PhasemaskTab);
            app.RotatephasemaskSlider.Limits = [0 1];
            app.RotatephasemaskSlider.ValueChangedFcn = createCallbackFcn(app, @RotatephasemaskSliderValueChanged, true);
            app.RotatephasemaskSlider.ValueChangingFcn = createCallbackFcn(app, @RotatephasemaskSliderValueChanging, true);
            app.RotatephasemaskSlider.BusyAction = 'cancel';
            app.RotatephasemaskSlider.Visible = 'off';
            app.RotatephasemaskSlider.Position = [136 195 150 3];

            % Create SectorSliderLabel
            app.SectorSliderLabel = uilabel(app.PhasemaskTab);
            app.SectorSliderLabel.Visible = 'off';
            app.SectorSliderLabel.Position = [16 125 40 22];
            app.SectorSliderLabel.Text = 'Sector';

            % Create SectorSlider
            app.SectorSlider = uislider(app.PhasemaskTab);
            app.SectorSlider.Limits = [0 360];
            app.SectorSlider.ValueChangedFcn = createCallbackFcn(app, @SectorSliderValueChanged, true);
            app.SectorSlider.ValueChangingFcn = createCallbackFcn(app, @SectorSliderValueChanging, true);
            app.SectorSlider.BusyAction = 'cancel';
            app.SectorSlider.Visible = 'off';
            app.SectorSlider.Position = [136 140 150 3];
            app.SectorSlider.Value = 180;

            % Create PhaseMaskOptionsLabel
            app.PhaseMaskOptionsLabel = uilabel(app.PhasemaskTab);
            app.PhaseMaskOptionsLabel.FontSize = 13;
            app.PhaseMaskOptionsLabel.FontWeight = 'bold';
            app.PhaseMaskOptionsLabel.Visible = 'off';
            app.PhaseMaskOptionsLabel.Position = [16 267 55 22];
            app.PhaseMaskOptionsLabel.Text = 'Options';

            % Create NumberfacetsSpinnerLabel
            app.NumberfacetsSpinnerLabel = uilabel(app.PhasemaskTab);
            app.NumberfacetsSpinnerLabel.Visible = 'off';
            app.NumberfacetsSpinnerLabel.Position = [16 124 92 22];
            app.NumberfacetsSpinnerLabel.Text = 'Number facets';

            % Create NumberfacetsSpinner
            app.NumberfacetsSpinner = uispinner(app.PhasemaskTab);
            app.NumberfacetsSpinner.ValueChangingFcn = createCallbackFcn(app, @NumberfacetsSpinnerValueChanging, true);
            app.NumberfacetsSpinner.Limits = [2 100];
            app.NumberfacetsSpinner.ValueDisplayFormat = '%.0f';
            app.NumberfacetsSpinner.ValueChangedFcn = createCallbackFcn(app, @NumberfacetsSpinnerValueChanged, true);
            app.NumberfacetsSpinner.BusyAction = 'cancel';
            app.NumberfacetsSpinner.Visible = 'off';
            app.NumberfacetsSpinner.Position = [136 123 56 22];
            app.NumberfacetsSpinner.Value = 2;

            % Create MaxShiftSliderLabel
            app.MaxShiftSliderLabel = uilabel(app.PhasemaskTab);
            app.MaxShiftSliderLabel.HorizontalAlignment = 'right';
            app.MaxShiftSliderLabel.Visible = 'off';
            app.MaxShiftSliderLabel.Position = [12 90 59 22];
            app.MaxShiftSliderLabel.Text = 'Max. Shift';

            % Create MaxShiftSlider
            app.MaxShiftSlider = uislider(app.PhasemaskTab);
            app.MaxShiftSlider.Limits = [0 6.28318530717959];
            app.MaxShiftSlider.MajorTicks = [0 3.14159265358979 6.28318530717959];
            app.MaxShiftSlider.MajorTickLabels = {'0', 'pi', '2pi'};
            app.MaxShiftSlider.ValueChangedFcn = createCallbackFcn(app, @MaxShiftSliderValueChanged, true);
            app.MaxShiftSlider.ValueChangingFcn = createCallbackFcn(app, @MaxShiftSliderValueChanging, true);
            app.MaxShiftSlider.Visible = 'off';
            app.MaxShiftSlider.Position = [136 100 101 3];
            app.MaxShiftSlider.Value = 6.28318530717959;

            % Create PhaseMaskShowplotCheckBox
            app.PhaseMaskShowplotCheckBox = uicheckbox(app.PhasemaskTab);
            app.PhaseMaskShowplotCheckBox.ValueChangedFcn = createCallbackFcn(app, @PhaseMaskShowplotCheckBoxValueChanged, true);
            app.PhaseMaskShowplotCheckBox.Text = 'Show plot';
            app.PhaseMaskShowplotCheckBox.Position = [120 352 76 22];

            % Create LoadCustomPhaseMaskButton
            app.LoadCustomPhaseMaskButton = uibutton(app.PhasemaskTab, 'push');
            app.LoadCustomPhaseMaskButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCustomPhaseMaskButtonPushed, true);
            app.LoadCustomPhaseMaskButton.Visible = 'off';
            app.LoadCustomPhaseMaskButton.Position = [301 315 114 23];
            app.LoadCustomPhaseMaskButton.Text = 'Load phase mask';

            % Create PhaseMaskFilepathLabel
            app.PhaseMaskFilepathLabel = uilabel(app.PhasemaskTab);
            app.PhaseMaskFilepathLabel.Visible = 'off';
            app.PhaseMaskFilepathLabel.Position = [136 283 352 22];
            app.PhaseMaskFilepathLabel.Text = '';

            % Create TransmissionTab
            app.TransmissionTab = uitab(app.TabGroup);
            app.TransmissionTab.Title = 'Transmission';

            % Create TransmissionLabel
            app.TransmissionLabel = uilabel(app.TransmissionTab);
            app.TransmissionLabel.FontSize = 13;
            app.TransmissionLabel.FontWeight = 'bold';
            app.TransmissionLabel.Position = [16 352 115 22];
            app.TransmissionLabel.Text = 'Transmission';

            % Create TransmissionMaskShowplotCheckBox
            app.TransmissionMaskShowplotCheckBox = uicheckbox(app.TransmissionTab);
            app.TransmissionMaskShowplotCheckBox.ValueChangedFcn = createCallbackFcn(app, @TransmissionMaskShowplotCheckBoxValueChanged, true);
            app.TransmissionMaskShowplotCheckBox.Text = 'Show plot';
            app.TransmissionMaskShowplotCheckBox.Position = [146 352 76 22];

            % Create TransmissionmaskDropDownLabel
            app.TransmissionmaskDropDownLabel = uilabel(app.TransmissionTab);
            app.TransmissionmaskDropDownLabel.Position = [16 315 108 22];
            app.TransmissionmaskDropDownLabel.Text = 'Transmission mask';

            % Create TransmissionmaskDropDown
            app.TransmissionmaskDropDown = uidropdown(app.TransmissionTab);
            app.TransmissionmaskDropDown.Items = {'None', 'Custom'};
            app.TransmissionmaskDropDown.ItemsData = {'none', 'Custom'};
            app.TransmissionmaskDropDown.ValueChangedFcn = createCallbackFcn(app, @TransmissionmaskDropDownValueChanged, true);
            app.TransmissionmaskDropDown.BusyAction = 'cancel';
            app.TransmissionmaskDropDown.Position = [146 315 145 22];
            app.TransmissionmaskDropDown.Value = 'none';

            % Create LoadCustomTransmissionMaskButton
            app.LoadCustomTransmissionMaskButton = uibutton(app.TransmissionTab, 'push');
            app.LoadCustomTransmissionMaskButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCustomTransmissionMaskButtonPushed, true);
            app.LoadCustomTransmissionMaskButton.Visible = 'off';
            app.LoadCustomTransmissionMaskButton.Position = [300 315 144 23];
            app.LoadCustomTransmissionMaskButton.Text = 'Load transmission mask';

            % Create TransmissionMaskFilepathLabel
            app.TransmissionMaskFilepathLabel = uilabel(app.TransmissionTab);
            app.TransmissionMaskFilepathLabel.Visible = 'off';
            app.TransmissionMaskFilepathLabel.Position = [146 283 352 22];
            app.TransmissionMaskFilepathLabel.Text = '';

            % Create OptionsTab
            app.OptionsTab = uitab(app.TabGroup);
            app.OptionsTab.Title = 'Options';
            app.OptionsTab.BusyAction = 'cancel';

            % Create SetcontrastButtonGroup
            app.SetcontrastButtonGroup = uibuttongroup(app.OptionsTab);
            app.SetcontrastButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @SetcontrastButtonGroupSelectionChanged, true);
            app.SetcontrastButtonGroup.BorderType = 'none';
            app.SetcontrastButtonGroup.BusyAction = 'cancel';
            app.SetcontrastButtonGroup.Position = [80 131 266 38];

            % Create IntensitylimitsButton
            app.IntensitylimitsButton = uiradiobutton(app.SetcontrastButtonGroup);
            app.IntensitylimitsButton.Tooltip = {'Colorbar limits set to 0 (lower limit) and maxium intensity (upper limit)'};
            app.IntensitylimitsButton.Text = 'Intensity limits';
            app.IntensitylimitsButton.Position = [10 11 127 22];
            app.IntensitylimitsButton.Value = true;

            % Create OptimizecontrastButton
            app.OptimizecontrastButton = uiradiobutton(app.SetcontrastButtonGroup);
            app.OptimizecontrastButton.BusyAction = 'cancel';
            app.OptimizecontrastButton.Tooltip = {'Colorbar limits set to minimum intensity (lower limit) and maximum intensity (upper value)'};
            app.OptimizecontrastButton.Text = 'Optimize contrast';
            app.OptimizecontrastButton.Position = [126 11 127 22];

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.OptionsTab);
            app.ColormapDropDownLabel.Position = [16 177 58 22];
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.OptionsTab);
            app.ColormapDropDown.Items = {'viridis', 'gray', 'parula', 'hot', 'jet', 'turbo'};
            app.ColormapDropDown.ValueChangedFcn = createCallbackFcn(app, @ColormapDropDownValueChanged, true);
            app.ColormapDropDown.BusyAction = 'cancel';
            app.ColormapDropDown.Position = [88 177 117 22];
            app.ColormapDropDown.Value = 'viridis';

            % Create ROIsidelengthEditFieldLabel
            app.ROIsidelengthEditFieldLabel = uilabel(app.OptionsTab);
            app.ROIsidelengthEditFieldLabel.BusyAction = 'cancel';
            app.ROIsidelengthEditFieldLabel.Tooltip = {''};
            app.ROIsidelengthEditFieldLabel.Position = [16 106 88 22];
            app.ROIsidelengthEditFieldLabel.Text = 'ROI side length';

            % Create ROIsidelengthEditField
            app.ROIsidelengthEditField = uieditfield(app.OptionsTab, 'numeric');
            app.ROIsidelengthEditField.Limits = [0 11];
            app.ROIsidelengthEditField.ValueDisplayFormat = '%11.4g µm';
            app.ROIsidelengthEditField.ValueChangedFcn = createCallbackFcn(app, @ROIsidelengthEditFieldValueChanged, true);
            app.ROIsidelengthEditField.BusyAction = 'cancel';
            app.ROIsidelengthEditField.Tooltip = {''};
            app.ROIsidelengthEditField.Position = [141 106 64 22];
            app.ROIsidelengthEditField.Value = 2.5;

            % Create PolarizedemissionchannelsCheckBox
            app.PolarizedemissionchannelsCheckBox = uicheckbox(app.OptionsTab);
            app.PolarizedemissionchannelsCheckBox.ValueChangedFcn = createCallbackFcn(app, @PolarizedemissionchannelsCheckBoxValueChanged, true);
            app.PolarizedemissionchannelsCheckBox.BusyAction = 'cancel';
            app.PolarizedemissionchannelsCheckBox.Tooltip = {'Note: Polarized emission images are shown without background or shot noise.'};
            app.PolarizedemissionchannelsCheckBox.Text = 'Polarized emission channels';
            app.PolarizedemissionchannelsCheckBox.Position = [16 247 199 22];

            % Create ShowPSFCheckBox
            app.ShowPSFCheckBox = uicheckbox(app.OptionsTab);
            app.ShowPSFCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowPSFCheckBoxValueChanged, true);
            app.ShowPSFCheckBox.BusyAction = 'cancel';
            app.ShowPSFCheckBox.Text = 'Point spread function';
            app.ShowPSFCheckBox.Position = [16 320 135 22];
            app.ShowPSFCheckBox.Value = true;

            % Create WindowselectionLabel
            app.WindowselectionLabel = uilabel(app.OptionsTab);
            app.WindowselectionLabel.FontSize = 13;
            app.WindowselectionLabel.FontWeight = 'bold';
            app.WindowselectionLabel.Tooltip = {'Select which windows are shown'};
            app.WindowselectionLabel.Position = [16 352 115 22];
            app.WindowselectionLabel.Text = 'Window selection';

            % Create PlotoptionsLabel
            app.PlotoptionsLabel = uilabel(app.OptionsTab);
            app.PlotoptionsLabel.FontSize = 13;
            app.PlotoptionsLabel.FontWeight = 'bold';
            app.PlotoptionsLabel.Tooltip = {'Select which windows are shown'};
            app.PlotoptionsLabel.Position = [16 208 80 22];
            app.PlotoptionsLabel.Text = 'Plot options';

            % Create ContrastLabel
            app.ContrastLabel = uilabel(app.OptionsTab);
            app.ContrastLabel.Position = [16 141 50 22];
            app.ContrastLabel.Text = 'Contrast';

            % Create PixelsperlateralaxisEditFieldLabel
            app.PixelsperlateralaxisEditFieldLabel = uilabel(app.OptionsTab);
            app.PixelsperlateralaxisEditFieldLabel.BusyAction = 'cancel';
            app.PixelsperlateralaxisEditFieldLabel.Tooltip = {''};
            app.PixelsperlateralaxisEditFieldLabel.Position = [16 76 118 22];
            app.PixelsperlateralaxisEditFieldLabel.Text = 'Pixels per lateral axis';

            % Create PixelsperlateralaxisEditField
            app.PixelsperlateralaxisEditField = uieditfield(app.OptionsTab, 'numeric');
            app.PixelsperlateralaxisEditField.Limits = [1 1000];
            app.PixelsperlateralaxisEditField.RoundFractionalValues = 'on';
            app.PixelsperlateralaxisEditField.ValueDisplayFormat = '%5d';
            app.PixelsperlateralaxisEditField.ValueChangedFcn = createCallbackFcn(app, @PixelsperlateralaxisEditFieldValueChanged, true);
            app.PixelsperlateralaxisEditField.BusyAction = 'cancel';
            app.PixelsperlateralaxisEditField.Tooltip = {''};
            app.PixelsperlateralaxisEditField.Position = [141 76 64 22];
            app.PixelsperlateralaxisEditField.Value = 25;

            % Create ShowPsf3DCheckBox
            app.ShowPsf3DCheckBox = uicheckbox(app.OptionsTab);
            app.ShowPsf3DCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowPsf3DCheckBoxValueChanged, true);
            app.ShowPsf3DCheckBox.BusyAction = 'cancel';
            app.ShowPsf3DCheckBox.Text = '3D';
            app.ShowPsf3DCheckBox.Position = [36 277 37 22];

            % Create ShowPsf2DCheckBox
            app.ShowPsf2DCheckBox = uicheckbox(app.OptionsTab);
            app.ShowPsf2DCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowPsf2DCheckBoxValueChanged, true);
            app.ShowPsf2DCheckBox.BusyAction = 'cancel';
            app.ShowPsf2DCheckBox.Text = '2D';
            app.ShowPsf2DCheckBox.Position = [36 298 37 22];
            app.ShowPsf2DCheckBox.Value = true;

            % Create NumberstepsEditFieldLabel
            app.NumberstepsEditFieldLabel = uilabel(app.OptionsTab);
            app.NumberstepsEditFieldLabel.BusyAction = 'cancel';
            app.NumberstepsEditFieldLabel.Enable = 'off';
            app.NumberstepsEditFieldLabel.Visible = 'off';
            app.NumberstepsEditFieldLabel.Tooltip = {''};
            app.NumberstepsEditFieldLabel.Position = [120 277 80 22];
            app.NumberstepsEditFieldLabel.Text = 'Number steps';

            % Create NumberstepsEditField
            app.NumberstepsEditField = uieditfield(app.OptionsTab, 'numeric');
            app.NumberstepsEditField.Limits = [1 1000];
            app.NumberstepsEditField.RoundFractionalValues = 'on';
            app.NumberstepsEditField.ValueDisplayFormat = '%5d';
            app.NumberstepsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberstepsEditFieldValueChanged, true);
            app.NumberstepsEditField.BusyAction = 'cancel';
            app.NumberstepsEditField.Enable = 'off';
            app.NumberstepsEditField.Visible = 'off';
            app.NumberstepsEditField.Tooltip = {''};
            app.NumberstepsEditField.Position = [210 277 36 22];
            app.NumberstepsEditField.Value = 25;

            % Create StepsizeEditFieldLabel
            app.StepsizeEditFieldLabel = uilabel(app.OptionsTab);
            app.StepsizeEditFieldLabel.BusyAction = 'cancel';
            app.StepsizeEditFieldLabel.Enable = 'off';
            app.StepsizeEditFieldLabel.Visible = 'off';
            app.StepsizeEditFieldLabel.Tooltip = {''};
            app.StepsizeEditFieldLabel.Position = [270 277 58 22];
            app.StepsizeEditFieldLabel.Text = 'Step size';

            % Create StepsizeEditField
            app.StepsizeEditField = uieditfield(app.OptionsTab, 'numeric');
            app.StepsizeEditField.Limits = [0 Inf];
            app.StepsizeEditField.ValueDisplayFormat = '%11.4g nm';
            app.StepsizeEditField.ValueChangedFcn = createCallbackFcn(app, @StepsizeEditFieldValueChanged, true);
            app.StepsizeEditField.BusyAction = 'cancel';
            app.StepsizeEditField.Enable = 'off';
            app.StepsizeEditField.Visible = 'off';
            app.StepsizeEditField.Position = [336 277 68 22];
            app.StepsizeEditField.Value = 100;

            % Create CalculateCramrRaoBoundCheckBox
            app.CalculateCramrRaoBoundCheckBox = uicheckbox(app.OptionsTab);
            app.CalculateCramrRaoBoundCheckBox.ValueChangedFcn = createCallbackFcn(app, @CalculateCramrRaoBoundCheckBoxValueChanged, true);
            app.CalculateCramrRaoBoundCheckBox.Text = 'Calculate Cramér Rao Bound';
            app.CalculateCramrRaoBoundCheckBox.Position = [388 175 179 22];

            % Create CRBOutputField
            app.CRBOutputField = uilabel(app.OptionsTab);
            app.CRBOutputField.Position = [407 122 137 47];
            app.CRBOutputField.Text = '';

            % Create CramrRaoBoundLabel
            app.CramrRaoBoundLabel = uilabel(app.OptionsTab);
            app.CramrRaoBoundLabel.FontSize = 13;
            app.CramrRaoBoundLabel.FontWeight = 'bold';
            app.CramrRaoBoundLabel.Tooltip = {'Select which windows are shown'};
            app.CramrRaoBoundLabel.Position = [388 208 125 22];
            app.CramrRaoBoundLabel.Text = 'Cramér-Rao Bound';

            % Create CalculatingLamp
            app.CalculatingLamp = uilamp(app.PSFsimulationUIFigure);
            app.CalculatingLamp.Interruptible = 'off';
            app.CalculatingLamp.Position = [602 435 20 20];

            % Show the figure after all components are created
            app.PSFsimulationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MainSimulationPSF_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PSFsimulationUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PSFsimulationUIFigure)
        end
    end
end