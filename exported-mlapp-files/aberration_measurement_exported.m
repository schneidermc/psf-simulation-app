classdef aberration_measurement_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PSFcharacterizationUIFigure     matlab.ui.Figure
        SetcolorlimitspersliceCheckBox  matlab.ui.control.CheckBox
        PloterrorButton                 matlab.ui.control.StateButton
        NumberZernikesEditField         matlab.ui.control.NumericEditField
        NumberZernikesEditFieldLabel    matlab.ui.control.Label
        parameterFile                   matlab.ui.control.EditField
        LoadparametersButton            matlab.ui.control.Button
        FitresultsLabel                 matlab.ui.control.Label
        PixelsizephysicalEditField      matlab.ui.control.NumericEditField
        PixelsizephysicalLabel          matlab.ui.control.Label
        RIimmersionEditField            matlab.ui.control.NumericEditField
        RIimmersionEditFieldLabel       matlab.ui.control.Label
        NAEditField                     matlab.ui.control.NumericEditField
        NAEditFieldLabel                matlab.ui.control.Label
        MagnificationEditField          matlab.ui.control.NumericEditField
        MagnificationEditFieldLabel     matlab.ui.control.Label
        RIsamplelayerEditField          matlab.ui.control.NumericEditField
        RIsamplelayerEditFieldLabel     matlab.ui.control.Label
        BeaddiameterEditField           matlab.ui.control.NumericEditField
        BeaddiameterEditFieldLabel      matlab.ui.control.Label
        EmissionwavelengthEditField     matlab.ui.control.NumericEditField
        EmissionwavelengthEditFieldLabel  matlab.ui.control.Label
        zincrementEditField             matlab.ui.control.NumericEditField
        zincrementEditFieldLabel        matlab.ui.control.Label
        IterationsEditField             matlab.ui.control.NumericEditField
        MaxiterationsLabel              matlab.ui.control.Label
        CalibrationsampleLabel          matlab.ui.control.Label
        MicroscopeparametersLabel       matlab.ui.control.Label
        text_stackfile                  matlab.ui.control.EditField
        ErrorLabel                      matlab.ui.control.Label
        Lamp                            matlab.ui.control.Lamp
        LoadzstackButton                matlab.ui.control.Button
        FitButton                       matlab.ui.control.Button
        zSliderSimulation               matlab.ui.control.Slider
        CalculatemodelButton            matlab.ui.control.Button
        zSliderMeasurement              matlab.ui.control.Slider
        FitZernikeaberrationsLabel      matlab.ui.control.Label
        SaveparametersButton            matlab.ui.control.Button
        SaveaberrationsButton           matlab.ui.control.Button
        SavetransmissionButton          matlab.ui.control.Button
        UIAxes_Transmission             matlab.ui.control.UIAxes
        UIAxes_Zernike                  matlab.ui.control.UIAxes
        UIAxes_PhaseAberration          matlab.ui.control.UIAxes
        UIAxesSimulationPsfXY           matlab.ui.control.UIAxes
        UIAxesSimulationPsfXZ           matlab.ui.control.UIAxes
        UIAxesMeasurementPsfXZ          matlab.ui.control.UIAxes
        UIAxesMeasurementPsfXY          matlab.ui.control.UIAxes
    end

    properties (Access = private)
        CallingApp

        obj = [] % objective information
        cam = [] % camera information
        zstack_folder = "./userData/zstack/"
        parameters_folder = "./userData/parameters/"
        aberrations_folder = "./userData/aberrations/"
        cam_loaded = false 
        obj_loaded = false 
        ux % spatial resolution
        
        stack % z-stack of experimental bead image
        stack_simu % z-stack of simulated bead image
        stack_error % error between experimental and simulated image
        PSF_xz % projection
        PSF_xy
        PSF_simu_xz
        PSF_simu_xy
        PSF_simu_xz_error
        PSF_simu_xy_error
        E_BFP % stack of 6 back focal plane fields
        dia_bead = 100e-9 % diameter of bead (in SI units [m])
        pupil % binary mask defining objective pupil in Fourier space
        dia_pupil = 64
        dz = 0.2e-6 % distance between two frames of the z-stack in µm
        m_focus = 9 % index of image where bead appears in focus
        lambda = 670e-9
        stack_loaded = 0 % boolean to check if stack has been loaded
        model_calculated = 0 % boolean to check if simu has been created
        info
        Nx = 50 % x-size of image stack
        Ny = 50
        Nz = 17 % depth of image stack (no. of images)
        RI_fluid = 1.33 % RI of bead solution
        Z_aberr % vector of retrieved zernike modes
        iter_max = 1000 % max. number of search iterations
        slice_norm = 1 % if set to 1, each z-slice is individually normalized (check if different exposure times have beeen used)
        Z_phase
        Z_amp
        modes % vector of Zernike modes used for fitting
        transmissionMask
    end
    
    methods (Access = private)
        function [colorMin, colorMax] = getColorLimitsPsf(app, slice)
            if nargin < 2
                slice = false;
            end
            if app.SetcolorlimitspersliceCheckBox.Value && slice
                sliceExp = app.stack(:,:,slice);
                if ~isempty(app.stack_simu)
                    sliceModel = app.stack_simu(:,:,slice);
                else
                    sliceModel = [];
                end
                colorMin = min([sliceExp(:); sliceModel(:)]);
                colorMax = max([sliceExp(:); sliceModel(:)]);
            else
                colorMin = min([app.stack_simu(:); app.stack(:)]);
                colorMax = max([app.stack_simu(:); app.stack(:)]);
            end
        end

        function [colorMin, colorMax] = getColorLimitsError(app, slice)
            if app.SetcolorlimitspersliceCheckBox.Value
                sliceValues = app.stack_error(:,:,slice);
                colorMin = min(sliceValues(:));
                colorMax = max(sliceValues(:));
            else
                colorMin = min(app.stack_error(:));
                colorMax = max(app.stack_error(:));
            end
        end

        function display_projection(app, axes, data_xz, Nx, Ny, Nz, ux, dz)
            x_axis = linspace(-(Nx-1)/2*ux*1e6, (Nx-1)/2*ux*1e6, Nx);
            z_axis = linspace(-(Nz-1)/2*dz*1e6, (Nz-1)/2*dz*1e6, Nz); % linspace(0, (Nz-1)*dz*1e6, Nz);%
            % Show xz projection
            h = axes; % select figure
            imagesc(h, x_axis, z_axis, data_xz);
            set(h,'DataAspectRatio',[dz, ux, ux])
            axis(h, "tight");
            colormap(h,app.CallingApp.ColormapDropDown.Value)
        end
               
        function display_single_image(app, axes, Nx, Ny, ux, data, cmin, cmax)
            h = axes; % select figure
            x_axis = linspace(-(Nx-1)/2*ux*1e6, (Nx-1)/2*ux*1e6, Nx);
            y_axis = linspace(-(Ny-1)/2*ux*1e6, (Ny-1)/2*ux*1e6, Ny);
            imagesc(h,x_axis, y_axis, data);
            set(h,'DataAspectRatio',[ux, ux, ux]);
            set(h,'YDir','normal')
            axis(h,"tight");
            colormap(h,app.CallingApp.ColormapDropDown.Value)
            caxis(h,[cmin,cmax])
        end

        function addColorbar(app, uiAxes)
            originalPosition = get(uiAxes, 'Position');
            cb = colorbar(uiAxes);
            set(uiAxes, 'Position', originalPosition);
            cbPosition = get(cb, 'Position');
            cbPosition(3) = 0.01; % Set the width
            set(cb, 'Position', cbPosition);
        end
        
        function drawZLineMeasurement(app,value)
            value = value+(1-app.Nz)/2*app.zincrementEditField.Value;
            persistent lineHandleMeasurement
            if isempty(lineHandleMeasurement) || ~isvalid(lineHandleMeasurement)
                lineHandleMeasurement = yline(app.UIAxesMeasurementPsfXZ, value,'w-');
            end
            % Update line position
            lineHandleMeasurement.Value = value; 
        end

        function drawZLineSimulation(app,value)
            value = value+(1-app.Nz)/2*app.zincrementEditField.Value;
            persistent lineHandleSimulation
            if isempty(lineHandleSimulation) || ~isvalid(lineHandleSimulation)
                lineHandleSimulation = yline(app.UIAxesSimulationPsfXZ, value,'w-');
            end
            % Update line position
            lineHandleSimulation.Value = value;
        end
        
        function display_result(~, axes, data)
            imagesc(axes, data);
            axis(axes, 'equal');
            axis(axes, 'tight');
            colorbar(axes);          
            axis(axes, "off");
        end

        function calculate_error(app)
            % Absolute error: abs(app.stack_simu - app.stack)
            % Relative error: abs(app.stack_simu - app.stack)./abs(app.stack_simu)
            app.stack_error = abs(app.stack_simu - app.stack);
            app.PSF_simu_xz_error = abs(app.PSF_simu_xz - app.PSF_xz);
            app.PSF_simu_xy_error = abs(app.PSF_simu_xy - app.PSF_xy);
        end
        
        function error = show_error(app, stack_A, stack_B)
            app.ErrorLabel.Visible = "on";
            error = std(stack_A(:) - stack_B(:));
            app.ErrorLabel.Text = "Error = " + num2str(error);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, mainapp)
            % Store main app object
            app.CallingApp = mainapp;
            app.EmissionwavelengthEditField.Value = app.CallingApp.EmissionWavelengthSpinner.Value;
            wavelength = Length(app.CallingApp.EmissionWavelengthSpinner.Value,'nm');
            app.lambda = wavelength.inMeter;
            
            app.RIsamplelayerEditField.Value = app.CallingApp.RefractiveIndexSpecimenSpinner.Value;
            app.RI_fluid = app.RIsamplelayerEditField.Value;
            
            % Objective
            app.MagnificationEditField.Value = app.CallingApp.MagnificationSpinner.Value;
            app.NAEditField.Value = app.CallingApp.ObjectiveNaSpinner.Value;
            app.RIimmersionEditField.Value = app.CallingApp.RefractiveIndexImmersionMediumSpinner.Value;
            app.obj = struct('name','Objective', ...
                'M',app.CallingApp.MagnificationSpinner.Value, ...
                'NA',app.CallingApp.ObjectiveNaSpinner.Value, ...
                'RI',app.CallingApp.RefractiveIndexImmersionMediumSpinner.Value);
            
            % Camera
            app.PixelsizephysicalEditField.Value = app.CallingApp.PixelsizePhysicalSpinner.Value;
            app.cam = struct('pixsize',app.CallingApp.PixelsizePhysicalSpinner.Value * 1e-6);

            app.UIAxesMeasurementPsfXY.PositionConstraint = 'innerposition';
            app.UIAxesMeasurementPsfXZ.PositionConstraint = 'innerposition';
            app.UIAxesSimulationPsfXY.PositionConstraint = 'innerposition';
            app.UIAxesSimulationPsfXZ.PositionConstraint = 'innerposition';

            app.modes = 2:app.NumberZernikesEditField.Value;
        end

        % Value changed function: BeaddiameterEditField
        function BeaddiameterEditFieldValueChanged(app, event)
            app.dia_bead = app.BeaddiameterEditField.Value * 1e-6;
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: zincrementEditField
        function zincrementEditFieldValueChanged(app, event)
            app.dz = app.zincrementEditField.Value*1e-6;
            if app.stack_loaded
                  display_projection(app, app.UIAxesMeasurementPsfXZ, app.PSF_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz)
                  flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
                  app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
            end
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: RIsamplelayerEditField
        function RIsamplelayerEditFieldValueChanged(app, event)
            app.RI_fluid = app.RIsamplelayerEditField.Value;
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: EmissionwavelengthEditField
        function EmissionwavelengthEditFieldValueChanged(app, event)
            app.lambda = app.EmissionwavelengthEditField.Value * 1e-9;
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: MagnificationEditField
        function MagnificationEditFieldValueChanged(app, event)
            app.MagnificationEditField.Value = event.Value;
            app.obj.M = app.MagnificationEditField.Value;
            app.ux = app.cam.pixsize/app.obj.M;
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: NAEditField
        function NAEditFieldValueChanged(app, event)
            app.NAEditField.Value = event.Value;
            app.obj.NA = app.NAEditField.Value;
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: RIimmersionEditField
        function RIimmersionEditFieldValueChanged(app, event)
            app.RIimmersionEditField.Value = event.Value;
            app.obj.RI = app.RIimmersionEditField.Value;
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Value changed function: PixelsizephysicalEditField
        function PixelsizephysicalEditFieldValueChanged(app, event)
            app.PixelsizephysicalEditField.Value = event.Value;
            app.cam.pixsize = app.PixelsizephysicalEditField.Value * 1e-6;
            app.ux = app.cam.pixsize/app.obj.M;
            if app.stack_loaded
                  display_projection(app, app.UIAxesMeasurementPsfXZ, app.PSF_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz)
                  flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
                  app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
            end
            if app.UIAxesSimulationPsfXZ.Visible == "on"
                CalculatemodelButtonPushed(app)
            end
        end

        % Button pushed function: LoadparametersButton
        function LoadparametersButtonValueChanged(app, event)
            [file, path] = uigetfile('*.mat');
            if file == 0
              % User clicked the Cancel button.
              return;
            end
            if isempty(who('-file', fullfile(path,file), "parameters"))
                errordlg('No parameters found in selected file!','Load error');
            else
                load(fullfile(path,file),"parameters");
                app.BeaddiameterEditField.Value = parameters.beadDiameter;
                app.zincrementEditField.Value = parameters.zIncrement;
                app.RIsamplelayerEditField.Value = parameters.RiFluid;
                app.EmissionwavelengthEditField.Value = parameters.emissionWavelength;
                app.MagnificationEditField.Value = parameters.magnification;
                app.NAEditField.Value = parameters.NA;
                app.RIimmersionEditField.Value = parameters.RiImmersion;
                app.PixelsizephysicalEditField.Value = parameters.pixelsizePhysical;
                app.parameterFile.Value = file;
                
                app.cam.pixsize = app.PixelsizephysicalEditField.Value * 1e-6;
                app.obj.M = app.MagnificationEditField.Value;
                app.ux = app.cam.pixsize/app.obj.M;
                app.obj.NA = app.NAEditField.Value;
                app.obj.RI = app.RIimmersionEditField.Value;
                if app.stack_loaded
                      display_projection(app, app.UIAxesMeasurementPsfXZ, app.PSF_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz)
                      flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
                      app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
                end
                if app.UIAxesSimulationPsfXZ.Visible == "on"
                    CalculatemodelButtonPushed(app)
                end
            end
        end

        % Button pushed function: LoadzstackButton
        function LoadzstackButtonPushed(app, event)
            % Loading stack of bead images
            [filename, pathname] = uigetfile('*.tif');
            if filename == 0
              % User clicked the Cancel button.
              return;
            end

            app.info = imfinfo([pathname, filename]);
            app.Nx = app.info(1).Height; %size of an image
            app.Ny = app.info(1).Width;
            app.Nz = length(app.info);
        
            app.text_stackfile.Value = filename; 

            app.stack = zeros(app.Nx, app.Ny, app.Nz); % initializing PSF stack
           
            for m = 1 : app.Nz
                % Read the kth image in this multipage tiff file.
                app.stack(:,:,m) = imread([pathname, filename], m);
            end	
             
            % Pre-processing steps
            [app.stack, app.m_focus] = centering(app.stack); % laterally centering
            app.stack = remove_BG(app.stack);
            %embed(app.stack, [app.Nx, app.Nx, app.Nz],0);
            
            app.ux = app.cam.pixsize/app.obj.M;

            % Normalization
            if app.slice_norm == 1
                app.stack = app.stack ./ sum(sum(app.stack,1),2); %normalizing each frame in the stack individually; this is to compensate for different exposure times
            end
            app.stack = app.stack/sum(app.stack(:)); 

            %app.PSF_xz = rot90(squeeze(max(app.stack,[],2))); %max. intensity projection
            app.PSF_xz = (squeeze(sum(app.stack,1)))'; %standard projection
            
            % Make items visible
            app.UIAxesMeasurementPsfXZ.Visible = "on";
            app.UIAxesMeasurementPsfXY.Visible = "on";
            %app.objectiveclosesttobeadLabel.Visible = "on";
            %app.objectivefarthestfrombeadLabel.Visible = "on";

            % Display the bead images
            display_projection(app, app.UIAxesMeasurementPsfXZ, app.PSF_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
            psfSlice = app.stack(:,:,round(app.zSliderMeasurement.Value));
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, psfSlice, min(psfSlice(:)), max(psfSlice(:)));
            app.addColorbar(app.UIAxesMeasurementPsfXY)
            
            app.stack_loaded = true; % flag that stack has been loaded
            app.CalculatemodelButton.Enable = "on";

            % Set slider parameters
            if app.zSliderSimulation.Visible == "off"
                app.zSliderMeasurement.Visible = "on";
                app.zSliderMeasurement.Limits = [1, app.Nz];
                app.zSliderMeasurement.MajorTicks = 1:app.Nz;
                app.zSliderMeasurement.MajorTickLabels = string(app.Nz + 1 - app.zSliderMeasurement.MajorTicks);
                app.zSliderMeasurement.MinorTicks = [];
                app.zSliderMeasurement.Visible = "on";
                app.zSliderMeasurement.Value = ceil(app.Nz/2);
                flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
                drawZLineMeasurement(app, (flippedValue-1)*app.zincrementEditField.Value);
            end
        end

        % Button pushed function: CalculatemodelButton
        function CalculatemodelButtonPushed(app, event)

            % Calculate simulated bead stack based on the experimental parameters
            F_bead = consider_bead_size(app.Nx, app.Nz, app.ux, app.dz, app.m_focus, app.dia_bead); % returns a low-pass-filter function to consider the finite bead-size
            [app.stack_simu, app.E_BFP] = simu_image_stack(app.dia_pupil, app.lambda, app.obj.NA, app.RI_fluid, app.obj.RI, app.dz, app.Nx, app.Nz, app.m_focus, app.ux, F_bead, app.slice_norm);
            app.PSF_simu_xz = (squeeze(sum(app.stack_simu,1)))'; % standard projection
            app.model_calculated = 1; 
            
            if app.stack_loaded == 1
                app.zSliderMeasurement.Enable = "off";
                app.zSliderMeasurement.Visible = "off";
            end

            % Display the simulated bead images
            if app.zSliderSimulation.Visible == "off"
                app.zSliderSimulation.Limits = [1, app.Nz];
                app.zSliderSimulation.MajorTicks = 1:app.Nz;
                app.zSliderSimulation.MajorTickLabels = string(app.Nz + 1 - app.zSliderSimulation.MajorTicks);
                app.zSliderSimulation.MinorTicks = [];
                app.zSliderSimulation.Visible = "on";
                
                % Make items visible
                app.zSliderSimulation.Visible = "on";
                app.UIAxesSimulationPsfXZ.Visible = "on";
                app.UIAxesSimulationPsfXY.Visible = "on";
            end

            if app.stack_loaded == 1
                val = round(app.zSliderMeasurement.Value);
                calculate_error(app)
                show_error(app, app.stack, app.stack_simu);
                app.FitButton.Enable = "on";
                app.PloterrorButton.Visible = "on";
                app.PloterrorButton.Enable = "on";
            else
                val = app.Nz;
            end
            flippedValue = app.Nz + 1 - val;
            app.zSliderSimulation.Value = val;
            app.PloterrorButton.Value = 0;

            display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
            psfSliceModel = app.stack_simu(:,:,val);
            psfSliceExperiment = app.stack(:,:,val);
            [colorMin, colorMax] = app.getColorLimitsPsf(val);
            display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, psfSliceModel, colorMin, colorMax);
            app.addColorbar(app.UIAxesSimulationPsfXY)
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, psfSliceExperiment, colorMin, colorMax);
            delete(app.UIAxesMeasurementPsfXY.Colorbar); % remove colorbar
            drawZLineSimulation(app, (flippedValue-1)*app.zincrementEditField.Value);
        end

        % Value changed function: IterationsEditField
        function IterationsEditFieldValueChanged(app, event)
            app.iter_max = app.IterationsEditField.Value;
        end

        % Button pushed function: FitButton
        function FitButtonPushed(app, event)
            app.Lamp.Visible = "on";
            
            z_vec=((1:app.Nz) - app.m_focus)*app.dz; %vector of z-positions 0=in-focus
            app.modes = 2:app.NumberZernikesEditField.Value; % vector of Zernike modes that are considered, 1 can be omitted (describes only a phase offset); Noll scheme is considered
            apo_no=4; %number of coefficients that are reserved for pupil apodization
                      %apo_no=1.... no apodization considered
                      %apo_no=2.... parabolic apodization considered
                      %apo_no=4.... up to 6th order apodization considered
                     
                                  
            %----pupil apodization modes: (they are following after the phase modes in "Zernike_stack"----
            [~,~,R_pupil, app.pupil]=create_coord(app.dia_pupil,1/app.dia_pupil*2,'exact');
            %---phase modes----: 
            Zernike_stack=flipud(ZernikeCalc(app.modes,ones(length(app.modes),1),app.pupil,'Noll'));
        
            Zernike_stack(:,:,length(app.modes)+1) = app.pupil;
            Zernike_stack(:,:,length(app.modes)+2) = -(R_pupil.^2) .* app.pupil;
            Zernike_stack(:,:,length(app.modes)+3) = -(R_pupil.^4) .* app.pupil;
            Zernike_stack(:,:,length(app.modes)+4) = -(R_pupil.^6) .* app.pupil; 

            Z_ini=[zeros(1,length(app.modes)) 1 zeros(1,apo_no-1)]; %initial estimates

            %----defining parameters for the fast chirped z-trafo----
            Nk = app.dia_pupil; %size of input field (in k-space)
            N_pad = app.Nx + Nk-0;
            k0 = 2*pi/app.lambda; 
            uk = 2*k0*app.obj.NA/app.dia_pupil;
            x=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5))';
            y=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5));
            ux_fft=2*pi/(Nk*uk); % fft resolution without padding
            r = app.ux/ux_fft; % required r-factor to meet the desired resolution at the given grid size N
            alpha = r*1/Nk;
            kernel = exp(-1i*alpha*pi*x.^2)*exp(-1i*alpha*pi*y.^2); % create quadratic convolution phase kernel; faster method
            F_kernel = fft2(ifftshift(kernel));

            [~, ~, Kr, ~] = create_coord(app.dia_pupil, uk, 'FFT');
            defocus = real(sqrt((app.obj.RI * k0)^2 - Kr.^2)); % defocus function 
            defocus = defocus - mean(defocus(app.pupil)); % subtract mean value (important for calculating inner products)
            
            F_bead = consider_bead_size(app.Nx, app.Nz, app.ux, app.dz, app.m_focus, app.dia_bead); %returns a low-pass-filter function to consider the finite bead-size
            algorithm = "trust-region-reflective"; % alternative: "levenberg-marquardt"

            

            %% FITTING
            resnorm_initial = sum(fun_errorcalc(Z_ini, app.stack, alpha, apo_no, kernel, F_kernel, defocus, z_vec, app.E_BFP, app.dia_pupil, Zernike_stack, F_bead, app.slice_norm).^2);

            d = uiprogressdlg(app.PSFcharacterizationUIFigure,'Title','Progress bar',...
                'Message','Fitting Zernike coefficients...','Cancelable','on');
            drawnow
            outputFunction = @(x,optimValues,state) updateProgress(x,optimValues,state,d,resnorm_initial, app.iter_max);
            options = optimoptions(@lsqnonlin,'Algorithm',algorithm ,'MaxFunctionEvaluations', app.iter_max,'OutputFcn',outputFunction);
            
            function stop = updateProgress(x, optimValues, state, d, resnorm_initial, maxIterations)
                % Update progress bar based on decrease in the residual norm and iteration number
            
                % Calculate progress based on residual norm
                resnorm_current = optimValues.resnorm;
                if resnorm_current < resnorm_initial
                    progress_resnorm = (resnorm_initial - resnorm_current) / resnorm_initial;
                else
                    progress_resnorm = 0;
                end
            
                % Calculate progress based on iterations
                progress_iteration = optimValues.funccount / maxIterations;
            
                % Combine the two progress measures, assuming equal weight
                progress = 0.3 * progress_resnorm + 0.7 * progress_iteration;
                progress = max(progress, 0.01);
                progress = min(progress, 1);
                d.Value = progress;
                stop = false; % This should only be true if you want to halt the algorithm.
            end

            app.Lamp.Color = "y";
            pause(0.01);
            LB = [-5 * ones(1,length(app.modes)), 0, zeros(1,apo_no-1)];
            UB = [+5 * ones(1,length(app.modes)), 2, ones(1,apo_no-1)];
            
            app.Z_aberr = lsqnonlin(@(Z_ini) fun_errorcalc(Z_ini, app.stack, alpha, apo_no, kernel, F_kernel, defocus, z_vec, app.E_BFP, app.dia_pupil, Zernike_stack, F_bead, app.slice_norm), Z_ini, LB, UB, options);
            close(d)
            app.Lamp.Color = "g";


            %% Separate outputs into phase and amplitude modes
            app.Z_phase = app.Z_aberr(1:end-apo_no); % coefficients for phase-aberrations
            app.Z_amp = app.Z_aberr(end-apo_no+1:end); % coefficients describing pupil apodization
            
            %phase = sum(ZernikeCalc(app.modes,[0 0 app.Z_phase(3:end)]',app.pupil,'Noll'),3); %calculating aberration phase function

            if apo_no==1 
                amp=1; 
            else
                tmp=zeros(1,1,length(app.Z_amp)); 
                tmp(1,1,:) = app.Z_amp;
                amp= sum(Zernike_stack(:,:,(end-apo_no+1):end).*repmat(tmp,[Nk,Nk,1]),3);
                amp=amp/max(amp(:));
            end

            % Remove spherical defocus
            defocus_coefs = ZernikeCalc(app.modes, defocus, app.pupil, 'Noll'); %zernike coefs of high-NA defocus
            a_def = dot(defocus_coefs/norm(defocus_coefs), app.Z_phase); 
            app.Z_phase = (app.Z_phase - a_def*defocus_coefs'/norm(defocus_coefs)); %defocus-free Zernike-coefs
            
            Z_tiptilt = app.Z_phase(1:2);
            app.Z_phase(1:2) = 0; %setting tip/tilt to zero
            phase = sum(flipud(ZernikeCalc(app.modes,app.Z_phase', app.pupil,'Noll')),3);
            app.Z_phase = app.Z_phase / (2*pi);
            
            %% Plot phase, transmission and Zernike coefficients
            % Make items visible
            app.FitresultsLabel.Visible = "on";
            app.UIAxes_PhaseAberration.Visible = "on";
            app.UIAxes_Transmission.Visible = "on";
            app.UIAxes_Zernike.Visible = "on";

            % Aberration
            aberrations = ZernikeAberrations(app.modes, app.Z_phase', app.dia_pupil);
            plot(aberrations, app.UIAxes_PhaseAberration);
            axis(app.UIAxes_PhaseAberration, 'equal');
            axis(app.UIAxes_PhaseAberration, 'tight');
            cb = colorbar(app.UIAxes_PhaseAberration);
            caxis(app.UIAxes_PhaseAberration, [-pi pi])
            cb.FontSize = 11;
            cb.Label.String = 'Phase shift';
            cb.Label.FontSize = 11;
            cb.Ticks = (-1:0.5:1)*pi;
            cb.TickLabels = {'-\pi','','0','','\pi'};
            set(app.UIAxes_PhaseAberration,'box','off')
            
            % Transmission
            fittedTransmissionMask = Transmission(amp.^2);
            plot(fittedTransmissionMask, app.UIAxes_Transmission);
            %imagesc(app.UIAxes_Transmission, fittedTransmissionMask.mask);
            axis(app.UIAxes_Transmission, 'equal');
            axis(app.UIAxes_Transmission, 'tight');
            cb = colorbar(app.UIAxes_Transmission);
            colormap(app.UIAxes_Transmission, gray)
            caxis(app.UIAxes_Transmission, [0 1])
            cb.FontSize = 11;
            cb.Label.String = 'Transmission';
            cb.Label.FontSize = 11;
            cb.Ticks = (0:0.25:1);
            cb.TickLabels = {'0','','0.5','','1'};
            set(app.UIAxes_Transmission,'box','off')

            app.transmissionMask = fittedTransmissionMask;
            
            % Zernikes
            stem(app.UIAxes_Zernike, app.modes(1:2), app.Z_phase(1:2), 'Color', 0.7*[1 1 1]); % tip/tilt (set to 0)
            hold(app.UIAxes_Zernike, 'on');
            stem(app.UIAxes_Zernike, app.modes(3:end), app.Z_phase(3:end), 'Color', [0, 0.4470, 0.7410]); % omit 2nd and 3rd (tip/tilt) as they are set to 0
            hold(app.UIAxes_Zernike, 'off');
            grid(app.UIAxes_Zernike, "on");
            
            % Check retrieved bead image
            [~, app.stack_simu,~]=fun_errorcalc(app.Z_aberr, app.stack/sum(app.stack(:)), alpha, apo_no, kernel, F_kernel, defocus, z_vec, app.E_BFP, app.dia_pupil, Zernike_stack, F_bead, app.slice_norm);
            
            % Display simulated bead images with aberrations taken into account
            app.PSF_simu_xz = (squeeze(sum(app.stack_simu,1)))'; %standard projection
            display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
            psfSliceModel = app.stack_simu(:,:,1);
            psfSliceExperiment = app.stack(:,:,1);
            [colorMin, colorMax] = app.getColorLimitsPsf(1);
            display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, psfSliceModel, colorMin, colorMax);
            app.addColorbar(app.UIAxesSimulationPsfXY)
            app.UIAxesSimulationPsfXZ.Title.String = 'Model';
            app.PloterrorButton.Text = 'Plot error';
            app.zSliderSimulation.Value = app.Nz;
            drawZLineMeasurement(app, 0);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, psfSliceExperiment, colorMin, colorMax);
            app.zSliderMeasurement.Value = app.Nz;
            drawZLineSimulation(app, 0);
            
            % Error
            if app.PloterrorButton.Value == 1
                app.PloterrorButton.Value = 0; % reset button
            end
            calculate_error(app)
            show_error(app, app.stack, app.stack_simu);

            app.SaveaberrationsButton.Enable = "on";
            app.SavetransmissionButton.Enable = "on";
        end

        % Value changed function: zSliderSimulation
        function zSliderSimulationValueChanged(app, event)
            app.zSliderSimulation.Value = round(event.Value);
            app.zSliderMeasurement.Value = app.zSliderSimulation.Value;
            flippedValue = app.Nz + 1 - app.zSliderSimulation.Value;
            
            [colorMin, colorMax] = app.getColorLimitsPsf(flippedValue);

            if app.PloterrorButton.Value
                % Show error plot
                [colorMinError, colorMaxError] = app.getColorLimitsError(flippedValue);
                data = app.stack_error(:,:,flippedValue);
                display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, data, colorMinError, colorMaxError)
                colormap(app.UIAxesSimulationPsfXZ,'hot')
                colormap(app.UIAxesSimulationPsfXY,'hot')
            else
                % Show fitted model
                data = app.stack_simu(:,:,flippedValue);
                display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, data, colorMin, colorMax)
                colormap(app.UIAxesSimulationPsfXZ,app.CallingApp.ColormapDropDown.Value)
                colormap(app.UIAxesSimulationPsfXY,app.CallingApp.ColormapDropDown.Value)
            end
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, app.stack(:,:,flippedValue), colorMin, colorMax)
            app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
            app.drawZLineSimulation((flippedValue-1)*app.zincrementEditField.Value)
        end

        % Value changing function: zSliderSimulation
        function zSliderSimulationValueChanging(app, event)
            zSliderSimulationValueChanged(app, event)
        end

        % Value changed function: zSliderMeasurement
        function zSliderMeasurementValueChanged(app, event)
            app.zSliderMeasurement.Value = round(event.Value);
            flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
            data = app.stack(:,:,flippedValue);
            [colorMin, colorMax] = app.getColorLimitsPsf(flippedValue);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, data,colorMin, colorMax)
            app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
        end

        % Value changing function: zSliderMeasurement
        function zSliderMeasurementValueChanging(app, event)
            app.zSliderMeasurement.Value = round(event.Value);
            flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
            data = app.stack(:,:,flippedValue);
            [colorMin, colorMax] = app.getColorLimitsPsf(flippedValue);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, data,colorMin, colorMax)
            app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
        end

        % Button pushed function: SaveparametersButton
        function SaveparametersButtonPushed(app, event)
            % Create struct of all parameters
            % Calibration sample
            parameters.beadDiameter = app.BeaddiameterEditField.Value; % µm
            parameters.zIncrement = app.zincrementEditField.Value; % µm
            parameters.RiFluid = app.RI_fluid;
            parameters.emissionWavelength = app.EmissionwavelengthEditField.Value;
            % Microscope parameters
            parameters.magnification = app.MagnificationEditField.Value;
            parameters.NA = app.NAEditField.Value;
            parameters.RiImmersion = app.RIimmersionEditField.Value;
            parameters.pixelsizePhysical = app.PixelsizephysicalEditField.Value;

            [file,path] = uiputfile('*.mat','Save parameter settings','parametersZernikeFit.mat');
            save(fullfile(path,file),"parameters")
        end

        % Button pushed function: SaveaberrationsButton
        function SaveaberrationsButtonPushed(app, event)
            % Create structure for saving aberration data
            % mat-file
            %aberrations.ZernikeIndices = app.modes;
            %aberrations.ZernikeCoefficients = app.Z_phase;
            %[file,path] = uiputfile('*.mat','Save results','aberrations.mat');
            %save(fullfile(path,file),"aberrations")
            % csv-file
            aberrations = [app.modes', app.Z_phase'];
            [file,path] = uiputfile('*.csv','Save results','aberrations.csv');
            writematrix(aberrations,fullfile(path,file),'Delimiter','tab')
        end

        % Button pushed function: SavetransmissionButton
        function SavetransmissionButtonPushed(app, event)
            transmission = app.transmissionMask.mask;            
            [file,path] = uiputfile('*.mat','Save results','transmission.mat');
            save(fullfile(path,file),"transmission")
        end

        % Value changed function: NumberZernikesEditField
        function NumberZernikesEditFieldValueChanged(app, event)
            app.NumberZernikesEditField.Value = event.Value;
            app.modes = 2:app.NumberZernikesEditField.Value;
        end

        % Value changed function: PloterrorButton
        function PloterrorButtonValueChanged(app, event)
            app.PloterrorButton.Value = event.Value;
            val = app.Nz + 1 - app.zSliderSimulation.Value;
            
            if app.PloterrorButton.Value
                app.PloterrorButton.Text = 'Plot model';
                % Show error plot
                [colorMinError, colorMaxError] = app.getColorLimitsError(val);
                display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz_error, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
                display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, app.stack_error(:,:,val), colorMinError, colorMaxError);
                app.UIAxesSimulationPsfXZ.Title.String = 'Model error';
                colormap(app.UIAxesSimulationPsfXZ,'hot')
                colormap(app.UIAxesSimulationPsfXY,'hot')
            else
                app.PloterrorButton.Text = 'Plot error';
                % Show fitted model
                [colorMin, colorMax] = app.getColorLimitsPsf(val);
                display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
                display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, app.stack_simu(:,:,val), colorMin, colorMax);
                app.UIAxesSimulationPsfXZ.Title.String = 'Model';
                colormap(app.UIAxesSimulationPsfXZ,app.CallingApp.ColormapDropDown.Value)
                colormap(app.UIAxesSimulationPsfXY,app.CallingApp.ColormapDropDown.Value)
            end
            app.drawZLineSimulation((val-1)*app.zincrementEditField.Value)
        end

        % Value changed function: SetcolorlimitspersliceCheckBox
        function SetcolorlimitspersliceCheckBoxValueChanged(app, event)
            app.SetcolorlimitspersliceCheckBox.Value = event.Value;
            
            if ~isempty(app.stack)
                % Update single slice plots
                if app.zSliderSimulation.Visible == "on"
                    val = app.Nz + 1 - app.zSliderSimulation.Value;
                else
                    val = app.Nz + 1 - app.zSliderMeasurement.Value;
                end
                
                [colorMin, colorMax] = app.getColorLimitsPsf(val);
    
                if app.PloterrorButton.Value
                    if ~isempty(app.stack_error)
                        % Show error plot
                        [colorMinError, colorMaxError] = app.getColorLimitsError(val);
                        display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz_error, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
                        display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, app.stack_error(:,:,val), colorMinError, colorMaxError);
                        app.UIAxesSimulationPsfXZ.Title.String = 'Model error';
                        colormap(app.UIAxesSimulationPsfXZ,'hot')
                        colormap(app.UIAxesSimulationPsfXY,'hot')
                        app.drawZLineSimulation((val-1)*app.zincrementEditField.Value)
                    end
                else
                    if ~isempty(app.stack_simu)
                        % Show fitted model
                        display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
                        display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, app.stack_simu(:,:,val), colorMin, colorMax);
                        app.UIAxesSimulationPsfXZ.Title.String = 'Model';
                        colormap(app.UIAxesSimulationPsfXZ,app.CallingApp.ColormapDropDown.Value)
                        colormap(app.UIAxesSimulationPsfXY,app.CallingApp.ColormapDropDown.Value)
                        app.drawZLineSimulation((val-1)*app.zincrementEditField.Value)
                    end
                end
                
                psfSliceExperiment = app.stack(:,:,val);
                display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, psfSliceExperiment, colorMin, colorMax);
                app.drawZLineMeasurement((val-1)*app.zincrementEditField.Value)
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PSFcharacterizationUIFigure and hide until all components are created
            app.PSFcharacterizationUIFigure = uifigure('Visible', 'off');
            app.PSFcharacterizationUIFigure.Position = [100 50 916 713];
            app.PSFcharacterizationUIFigure.Name = 'PSF characterization';

            % Create UIAxesMeasurementPsfXY
            app.UIAxesMeasurementPsfXY = uiaxes(app.PSFcharacterizationUIFigure);
            xlabel(app.UIAxesMeasurementPsfXY, 'X / µm')
            ylabel(app.UIAxesMeasurementPsfXY, 'Y / µm')
            zlabel(app.UIAxesMeasurementPsfXY, 'Z')
            app.UIAxesMeasurementPsfXY.FontSize = 10;
            app.UIAxesMeasurementPsfXY.Visible = 'off';
            app.UIAxesMeasurementPsfXY.Position = [37 15 200 200];

            % Create UIAxesMeasurementPsfXZ
            app.UIAxesMeasurementPsfXZ = uiaxes(app.PSFcharacterizationUIFigure);
            title(app.UIAxesMeasurementPsfXZ, 'Measurement')
            ylabel(app.UIAxesMeasurementPsfXZ, 'Z / µm')
            zlabel(app.UIAxesMeasurementPsfXZ, 'Z')
            app.UIAxesMeasurementPsfXZ.FontSize = 10;
            app.UIAxesMeasurementPsfXZ.Visible = 'off';
            app.UIAxesMeasurementPsfXZ.Position = [37 222 200 270];

            % Create UIAxesSimulationPsfXZ
            app.UIAxesSimulationPsfXZ = uiaxes(app.PSFcharacterizationUIFigure);
            title(app.UIAxesSimulationPsfXZ, 'Simulation')
            ylabel(app.UIAxesSimulationPsfXZ, 'Z / µm')
            zlabel(app.UIAxesSimulationPsfXZ, 'Z')
            app.UIAxesSimulationPsfXZ.FontSize = 10;
            app.UIAxesSimulationPsfXZ.Visible = 'off';
            app.UIAxesSimulationPsfXZ.Position = [260 222 200 270];

            % Create UIAxesSimulationPsfXY
            app.UIAxesSimulationPsfXY = uiaxes(app.PSFcharacterizationUIFigure);
            xlabel(app.UIAxesSimulationPsfXY, 'X / µm')
            ylabel(app.UIAxesSimulationPsfXY, 'Y / µm')
            zlabel(app.UIAxesSimulationPsfXY, 'Z')
            app.UIAxesSimulationPsfXY.FontSize = 10;
            app.UIAxesSimulationPsfXY.Visible = 'off';
            app.UIAxesSimulationPsfXY.Position = [260 15 200 200];

            % Create UIAxes_PhaseAberration
            app.UIAxes_PhaseAberration = uiaxes(app.PSFcharacterizationUIFigure);
            title(app.UIAxes_PhaseAberration, 'Aberrations')
            app.UIAxes_PhaseAberration.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_PhaseAberration.XTick = [];
            app.UIAxes_PhaseAberration.YTick = [];
            app.UIAxes_PhaseAberration.ZTick = [];
            app.UIAxes_PhaseAberration.Box = 'on';
            app.UIAxes_PhaseAberration.Visible = 'off';
            app.UIAxes_PhaseAberration.Position = [549 289 160 160];

            % Create UIAxes_Zernike
            app.UIAxes_Zernike = uiaxes(app.PSFcharacterizationUIFigure);
            title(app.UIAxes_Zernike, 'Zernike modes')
            xlabel(app.UIAxes_Zernike, 'Noll index')
            ylabel(app.UIAxes_Zernike, 'RMS wavefront error [λ]')
            zlabel(app.UIAxes_Zernike, 'Z')
            app.UIAxes_Zernike.YGrid = 'on';
            app.UIAxes_Zernike.Box = 'on';
            app.UIAxes_Zernike.Visible = 'off';
            app.UIAxes_Zernike.Position = [549 96 345 185];

            % Create UIAxes_Transmission
            app.UIAxes_Transmission = uiaxes(app.PSFcharacterizationUIFigure);
            title(app.UIAxes_Transmission, 'Transmission')
            app.UIAxes_Transmission.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_Transmission.Colormap = [0 0 0;0.00392156862745098 0.00392156862745098 0.00392156862745098;0.00784313725490196 0.00784313725490196 0.00784313725490196;0.0117647058823529 0.0117647058823529 0.0117647058823529;0.0156862745098039 0.0156862745098039 0.0156862745098039;0.0196078431372549 0.0196078431372549 0.0196078431372549;0.0235294117647059 0.0235294117647059 0.0235294117647059;0.0274509803921569 0.0274509803921569 0.0274509803921569;0.0313725490196078 0.0313725490196078 0.0313725490196078;0.0352941176470588 0.0352941176470588 0.0352941176470588;0.0392156862745098 0.0392156862745098 0.0392156862745098;0.0431372549019608 0.0431372549019608 0.0431372549019608;0.0470588235294118 0.0470588235294118 0.0470588235294118;0.0509803921568627 0.0509803921568627 0.0509803921568627;0.0549019607843137 0.0549019607843137 0.0549019607843137;0.0588235294117647 0.0588235294117647 0.0588235294117647;0.0627450980392157 0.0627450980392157 0.0627450980392157;0.0666666666666667 0.0666666666666667 0.0666666666666667;0.0705882352941176 0.0705882352941176 0.0705882352941176;0.0745098039215686 0.0745098039215686 0.0745098039215686;0.0784313725490196 0.0784313725490196 0.0784313725490196;0.0823529411764706 0.0823529411764706 0.0823529411764706;0.0862745098039216 0.0862745098039216 0.0862745098039216;0.0901960784313725 0.0901960784313725 0.0901960784313725;0.0941176470588235 0.0941176470588235 0.0941176470588235;0.0980392156862745 0.0980392156862745 0.0980392156862745;0.101960784313725 0.101960784313725 0.101960784313725;0.105882352941176 0.105882352941176 0.105882352941176;0.109803921568627 0.109803921568627 0.109803921568627;0.113725490196078 0.113725490196078 0.113725490196078;0.117647058823529 0.117647058823529 0.117647058823529;0.12156862745098 0.12156862745098 0.12156862745098;0.125490196078431 0.125490196078431 0.125490196078431;0.129411764705882 0.129411764705882 0.129411764705882;0.133333333333333 0.133333333333333 0.133333333333333;0.137254901960784 0.137254901960784 0.137254901960784;0.141176470588235 0.141176470588235 0.141176470588235;0.145098039215686 0.145098039215686 0.145098039215686;0.149019607843137 0.149019607843137 0.149019607843137;0.152941176470588 0.152941176470588 0.152941176470588;0.156862745098039 0.156862745098039 0.156862745098039;0.16078431372549 0.16078431372549 0.16078431372549;0.164705882352941 0.164705882352941 0.164705882352941;0.168627450980392 0.168627450980392 0.168627450980392;0.172549019607843 0.172549019607843 0.172549019607843;0.176470588235294 0.176470588235294 0.176470588235294;0.180392156862745 0.180392156862745 0.180392156862745;0.184313725490196 0.184313725490196 0.184313725490196;0.188235294117647 0.188235294117647 0.188235294117647;0.192156862745098 0.192156862745098 0.192156862745098;0.196078431372549 0.196078431372549 0.196078431372549;0.2 0.2 0.2;0.203921568627451 0.203921568627451 0.203921568627451;0.207843137254902 0.207843137254902 0.207843137254902;0.211764705882353 0.211764705882353 0.211764705882353;0.215686274509804 0.215686274509804 0.215686274509804;0.219607843137255 0.219607843137255 0.219607843137255;0.223529411764706 0.223529411764706 0.223529411764706;0.227450980392157 0.227450980392157 0.227450980392157;0.231372549019608 0.231372549019608 0.231372549019608;0.235294117647059 0.235294117647059 0.235294117647059;0.23921568627451 0.23921568627451 0.23921568627451;0.243137254901961 0.243137254901961 0.243137254901961;0.247058823529412 0.247058823529412 0.247058823529412;0.250980392156863 0.250980392156863 0.250980392156863;0.254901960784314 0.254901960784314 0.254901960784314;0.258823529411765 0.258823529411765 0.258823529411765;0.262745098039216 0.262745098039216 0.262745098039216;0.266666666666667 0.266666666666667 0.266666666666667;0.270588235294118 0.270588235294118 0.270588235294118;0.274509803921569 0.274509803921569 0.274509803921569;0.27843137254902 0.27843137254902 0.27843137254902;0.282352941176471 0.282352941176471 0.282352941176471;0.286274509803922 0.286274509803922 0.286274509803922;0.290196078431373 0.290196078431373 0.290196078431373;0.294117647058824 0.294117647058824 0.294117647058824;0.298039215686275 0.298039215686275 0.298039215686275;0.301960784313725 0.301960784313725 0.301960784313725;0.305882352941176 0.305882352941176 0.305882352941176;0.309803921568627 0.309803921568627 0.309803921568627;0.313725490196078 0.313725490196078 0.313725490196078;0.317647058823529 0.317647058823529 0.317647058823529;0.32156862745098 0.32156862745098 0.32156862745098;0.325490196078431 0.325490196078431 0.325490196078431;0.329411764705882 0.329411764705882 0.329411764705882;0.333333333333333 0.333333333333333 0.333333333333333;0.337254901960784 0.337254901960784 0.337254901960784;0.341176470588235 0.341176470588235 0.341176470588235;0.345098039215686 0.345098039215686 0.345098039215686;0.349019607843137 0.349019607843137 0.349019607843137;0.352941176470588 0.352941176470588 0.352941176470588;0.356862745098039 0.356862745098039 0.356862745098039;0.36078431372549 0.36078431372549 0.36078431372549;0.364705882352941 0.364705882352941 0.364705882352941;0.368627450980392 0.368627450980392 0.368627450980392;0.372549019607843 0.372549019607843 0.372549019607843;0.376470588235294 0.376470588235294 0.376470588235294;0.380392156862745 0.380392156862745 0.380392156862745;0.384313725490196 0.384313725490196 0.384313725490196;0.388235294117647 0.388235294117647 0.388235294117647;0.392156862745098 0.392156862745098 0.392156862745098;0.396078431372549 0.396078431372549 0.396078431372549;0.4 0.4 0.4;0.403921568627451 0.403921568627451 0.403921568627451;0.407843137254902 0.407843137254902 0.407843137254902;0.411764705882353 0.411764705882353 0.411764705882353;0.415686274509804 0.415686274509804 0.415686274509804;0.419607843137255 0.419607843137255 0.419607843137255;0.423529411764706 0.423529411764706 0.423529411764706;0.427450980392157 0.427450980392157 0.427450980392157;0.431372549019608 0.431372549019608 0.431372549019608;0.435294117647059 0.435294117647059 0.435294117647059;0.43921568627451 0.43921568627451 0.43921568627451;0.443137254901961 0.443137254901961 0.443137254901961;0.447058823529412 0.447058823529412 0.447058823529412;0.450980392156863 0.450980392156863 0.450980392156863;0.454901960784314 0.454901960784314 0.454901960784314;0.458823529411765 0.458823529411765 0.458823529411765;0.462745098039216 0.462745098039216 0.462745098039216;0.466666666666667 0.466666666666667 0.466666666666667;0.470588235294118 0.470588235294118 0.470588235294118;0.474509803921569 0.474509803921569 0.474509803921569;0.47843137254902 0.47843137254902 0.47843137254902;0.482352941176471 0.482352941176471 0.482352941176471;0.486274509803922 0.486274509803922 0.486274509803922;0.490196078431373 0.490196078431373 0.490196078431373;0.494117647058824 0.494117647058824 0.494117647058824;0.498039215686275 0.498039215686275 0.498039215686275;0.501960784313725 0.501960784313725 0.501960784313725;0.505882352941176 0.505882352941176 0.505882352941176;0.509803921568627 0.509803921568627 0.509803921568627;0.513725490196078 0.513725490196078 0.513725490196078;0.517647058823529 0.517647058823529 0.517647058823529;0.52156862745098 0.52156862745098 0.52156862745098;0.525490196078431 0.525490196078431 0.525490196078431;0.529411764705882 0.529411764705882 0.529411764705882;0.533333333333333 0.533333333333333 0.533333333333333;0.537254901960784 0.537254901960784 0.537254901960784;0.541176470588235 0.541176470588235 0.541176470588235;0.545098039215686 0.545098039215686 0.545098039215686;0.549019607843137 0.549019607843137 0.549019607843137;0.552941176470588 0.552941176470588 0.552941176470588;0.556862745098039 0.556862745098039 0.556862745098039;0.56078431372549 0.56078431372549 0.56078431372549;0.564705882352941 0.564705882352941 0.564705882352941;0.568627450980392 0.568627450980392 0.568627450980392;0.572549019607843 0.572549019607843 0.572549019607843;0.576470588235294 0.576470588235294 0.576470588235294;0.580392156862745 0.580392156862745 0.580392156862745;0.584313725490196 0.584313725490196 0.584313725490196;0.588235294117647 0.588235294117647 0.588235294117647;0.592156862745098 0.592156862745098 0.592156862745098;0.596078431372549 0.596078431372549 0.596078431372549;0.6 0.6 0.6;0.603921568627451 0.603921568627451 0.603921568627451;0.607843137254902 0.607843137254902 0.607843137254902;0.611764705882353 0.611764705882353 0.611764705882353;0.615686274509804 0.615686274509804 0.615686274509804;0.619607843137255 0.619607843137255 0.619607843137255;0.623529411764706 0.623529411764706 0.623529411764706;0.627450980392157 0.627450980392157 0.627450980392157;0.631372549019608 0.631372549019608 0.631372549019608;0.635294117647059 0.635294117647059 0.635294117647059;0.63921568627451 0.63921568627451 0.63921568627451;0.643137254901961 0.643137254901961 0.643137254901961;0.647058823529412 0.647058823529412 0.647058823529412;0.650980392156863 0.650980392156863 0.650980392156863;0.654901960784314 0.654901960784314 0.654901960784314;0.658823529411765 0.658823529411765 0.658823529411765;0.662745098039216 0.662745098039216 0.662745098039216;0.666666666666667 0.666666666666667 0.666666666666667;0.670588235294118 0.670588235294118 0.670588235294118;0.674509803921569 0.674509803921569 0.674509803921569;0.67843137254902 0.67843137254902 0.67843137254902;0.682352941176471 0.682352941176471 0.682352941176471;0.686274509803922 0.686274509803922 0.686274509803922;0.690196078431373 0.690196078431373 0.690196078431373;0.694117647058824 0.694117647058824 0.694117647058824;0.698039215686274 0.698039215686274 0.698039215686274;0.701960784313725 0.701960784313725 0.701960784313725;0.705882352941177 0.705882352941177 0.705882352941177;0.709803921568627 0.709803921568627 0.709803921568627;0.713725490196078 0.713725490196078 0.713725490196078;0.717647058823529 0.717647058823529 0.717647058823529;0.72156862745098 0.72156862745098 0.72156862745098;0.725490196078431 0.725490196078431 0.725490196078431;0.729411764705882 0.729411764705882 0.729411764705882;0.733333333333333 0.733333333333333 0.733333333333333;0.737254901960784 0.737254901960784 0.737254901960784;0.741176470588235 0.741176470588235 0.741176470588235;0.745098039215686 0.745098039215686 0.745098039215686;0.749019607843137 0.749019607843137 0.749019607843137;0.752941176470588 0.752941176470588 0.752941176470588;0.756862745098039 0.756862745098039 0.756862745098039;0.76078431372549 0.76078431372549 0.76078431372549;0.764705882352941 0.764705882352941 0.764705882352941;0.768627450980392 0.768627450980392 0.768627450980392;0.772549019607843 0.772549019607843 0.772549019607843;0.776470588235294 0.776470588235294 0.776470588235294;0.780392156862745 0.780392156862745 0.780392156862745;0.784313725490196 0.784313725490196 0.784313725490196;0.788235294117647 0.788235294117647 0.788235294117647;0.792156862745098 0.792156862745098 0.792156862745098;0.796078431372549 0.796078431372549 0.796078431372549;0.8 0.8 0.8;0.803921568627451 0.803921568627451 0.803921568627451;0.807843137254902 0.807843137254902 0.807843137254902;0.811764705882353 0.811764705882353 0.811764705882353;0.815686274509804 0.815686274509804 0.815686274509804;0.819607843137255 0.819607843137255 0.819607843137255;0.823529411764706 0.823529411764706 0.823529411764706;0.827450980392157 0.827450980392157 0.827450980392157;0.831372549019608 0.831372549019608 0.831372549019608;0.835294117647059 0.835294117647059 0.835294117647059;0.83921568627451 0.83921568627451 0.83921568627451;0.843137254901961 0.843137254901961 0.843137254901961;0.847058823529412 0.847058823529412 0.847058823529412;0.850980392156863 0.850980392156863 0.850980392156863;0.854901960784314 0.854901960784314 0.854901960784314;0.858823529411765 0.858823529411765 0.858823529411765;0.862745098039216 0.862745098039216 0.862745098039216;0.866666666666667 0.866666666666667 0.866666666666667;0.870588235294118 0.870588235294118 0.870588235294118;0.874509803921569 0.874509803921569 0.874509803921569;0.87843137254902 0.87843137254902 0.87843137254902;0.882352941176471 0.882352941176471 0.882352941176471;0.886274509803922 0.886274509803922 0.886274509803922;0.890196078431372 0.890196078431372 0.890196078431372;0.894117647058824 0.894117647058824 0.894117647058824;0.898039215686275 0.898039215686275 0.898039215686275;0.901960784313726 0.901960784313726 0.901960784313726;0.905882352941176 0.905882352941176 0.905882352941176;0.909803921568627 0.909803921568627 0.909803921568627;0.913725490196078 0.913725490196078 0.913725490196078;0.917647058823529 0.917647058823529 0.917647058823529;0.92156862745098 0.92156862745098 0.92156862745098;0.925490196078431 0.925490196078431 0.925490196078431;0.929411764705882 0.929411764705882 0.929411764705882;0.933333333333333 0.933333333333333 0.933333333333333;0.937254901960784 0.937254901960784 0.937254901960784;0.941176470588235 0.941176470588235 0.941176470588235;0.945098039215686 0.945098039215686 0.945098039215686;0.949019607843137 0.949019607843137 0.949019607843137;0.952941176470588 0.952941176470588 0.952941176470588;0.956862745098039 0.956862745098039 0.956862745098039;0.96078431372549 0.96078431372549 0.96078431372549;0.964705882352941 0.964705882352941 0.964705882352941;0.968627450980392 0.968627450980392 0.968627450980392;0.972549019607843 0.972549019607843 0.972549019607843;0.976470588235294 0.976470588235294 0.976470588235294;0.980392156862745 0.980392156862745 0.980392156862745;0.984313725490196 0.984313725490196 0.984313725490196;0.988235294117647 0.988235294117647 0.988235294117647;0.992156862745098 0.992156862745098 0.992156862745098;0.996078431372549 0.996078431372549 0.996078431372549;1 1 1];
            app.UIAxes_Transmission.XTick = [];
            app.UIAxes_Transmission.YTick = [];
            app.UIAxes_Transmission.ZTick = [];
            app.UIAxes_Transmission.Box = 'on';
            app.UIAxes_Transmission.Visible = 'off';
            app.UIAxes_Transmission.Position = [729 289 165 160];

            % Create SavetransmissionButton
            app.SavetransmissionButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.SavetransmissionButton.ButtonPushedFcn = createCallbackFcn(app, @SavetransmissionButtonPushed, true);
            app.SavetransmissionButton.Enable = 'off';
            app.SavetransmissionButton.Position = [781 15 113 23];
            app.SavetransmissionButton.Text = 'Save transmission';

            % Create SaveaberrationsButton
            app.SaveaberrationsButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.SaveaberrationsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveaberrationsButtonPushed, true);
            app.SaveaberrationsButton.Enable = 'off';
            app.SaveaberrationsButton.Position = [665 15 106 23];
            app.SaveaberrationsButton.Text = 'Save aberrations';

            % Create SaveparametersButton
            app.SaveparametersButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.SaveparametersButton.ButtonPushedFcn = createCallbackFcn(app, @SaveparametersButtonPushed, true);
            app.SaveparametersButton.Position = [548 15 106 23];
            app.SaveparametersButton.Text = 'Save parameters';

            % Create FitZernikeaberrationsLabel
            app.FitZernikeaberrationsLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.FitZernikeaberrationsLabel.FontSize = 14;
            app.FitZernikeaberrationsLabel.FontWeight = 'bold';
            app.FitZernikeaberrationsLabel.Position = [18 682 156 22];
            app.FitZernikeaberrationsLabel.Text = 'Fit Zernike aberrations';

            % Create zSliderMeasurement
            app.zSliderMeasurement = uislider(app.PSFcharacterizationUIFigure);
            app.zSliderMeasurement.Orientation = 'vertical';
            app.zSliderMeasurement.ValueChangedFcn = createCallbackFcn(app, @zSliderMeasurementValueChanged, true);
            app.zSliderMeasurement.ValueChangingFcn = createCallbackFcn(app, @zSliderMeasurementValueChanging, true);
            app.zSliderMeasurement.BusyAction = 'cancel';
            app.zSliderMeasurement.Visible = 'off';
            app.zSliderMeasurement.FontSize = 10;
            app.zSliderMeasurement.Position = [240 284 3 150];
            app.zSliderMeasurement.Value = 1;

            % Create CalculatemodelButton
            app.CalculatemodelButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.CalculatemodelButton.ButtonPushedFcn = createCallbackFcn(app, @CalculatemodelButtonPushed, true);
            app.CalculatemodelButton.VerticalAlignment = 'top';
            app.CalculatemodelButton.Enable = 'off';
            app.CalculatemodelButton.Position = [548 568 110 23];
            app.CalculatemodelButton.Text = 'Calculate model';

            % Create zSliderSimulation
            app.zSliderSimulation = uislider(app.PSFcharacterizationUIFigure);
            app.zSliderSimulation.Orientation = 'vertical';
            app.zSliderSimulation.ValueChangedFcn = createCallbackFcn(app, @zSliderSimulationValueChanged, true);
            app.zSliderSimulation.ValueChangingFcn = createCallbackFcn(app, @zSliderSimulationValueChanging, true);
            app.zSliderSimulation.BusyAction = 'cancel';
            app.zSliderSimulation.Visible = 'off';
            app.zSliderSimulation.FontSize = 10;
            app.zSliderSimulation.Position = [468 284 3 150];

            % Create FitButton
            app.FitButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.FitButton.ButtonPushedFcn = createCallbackFcn(app, @FitButtonPushed, true);
            app.FitButton.BackgroundColor = [0.9608 0.9608 0.9608];
            app.FitButton.Enable = 'off';
            app.FitButton.Position = [830 539 39 23];
            app.FitButton.Text = 'Fit';

            % Create LoadzstackButton
            app.LoadzstackButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.LoadzstackButton.ButtonPushedFcn = createCallbackFcn(app, @LoadzstackButtonPushed, true);
            app.LoadzstackButton.Position = [549 596 110 23];
            app.LoadzstackButton.Text = 'Load z-stack';

            % Create Lamp
            app.Lamp = uilamp(app.PSFcharacterizationUIFigure);
            app.Lamp.Visible = 'off';
            app.Lamp.Position = [874 541 20 20];
            app.Lamp.Color = [0.651 0.651 0.651];

            % Create ErrorLabel
            app.ErrorLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.ErrorLabel.Visible = 'off';
            app.ErrorLabel.Position = [667 568 148 22];
            app.ErrorLabel.Text = 'Error = ';

            % Create text_stackfile
            app.text_stackfile = uieditfield(app.PSFcharacterizationUIFigure, 'text');
            app.text_stackfile.Editable = 'off';
            app.text_stackfile.Position = [667 597 227 22];

            % Create MicroscopeparametersLabel
            app.MicroscopeparametersLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.MicroscopeparametersLabel.FontWeight = 'bold';
            app.MicroscopeparametersLabel.Position = [248 653 140 22];
            app.MicroscopeparametersLabel.Text = 'Microscope parameters';

            % Create CalibrationsampleLabel
            app.CalibrationsampleLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.CalibrationsampleLabel.FontWeight = 'bold';
            app.CalibrationsampleLabel.Position = [18 653 112 22];
            app.CalibrationsampleLabel.Text = 'Calibration sample';

            % Create MaxiterationsLabel
            app.MaxiterationsLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.MaxiterationsLabel.HorizontalAlignment = 'right';
            app.MaxiterationsLabel.Position = [710 539 58 22];
            app.MaxiterationsLabel.Text = 'Iterations ';

            % Create IterationsEditField
            app.IterationsEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.IterationsEditField.ValueChangedFcn = createCallbackFcn(app, @IterationsEditFieldValueChanged, true);
            app.IterationsEditField.Position = [772 539 43 22];
            app.IterationsEditField.Value = 1000;

            % Create zincrementEditFieldLabel
            app.zincrementEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.zincrementEditFieldLabel.HorizontalAlignment = 'right';
            app.zincrementEditFieldLabel.Position = [68 597 68 22];
            app.zincrementEditFieldLabel.Text = 'z-increment';

            % Create zincrementEditField
            app.zincrementEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.zincrementEditField.Limits = [0.001 Inf];
            app.zincrementEditField.ValueDisplayFormat = '%11.4gµm';
            app.zincrementEditField.ValueChangedFcn = createCallbackFcn(app, @zincrementEditFieldValueChanged, true);
            app.zincrementEditField.Position = [145 597 58 22];
            app.zincrementEditField.Value = 0.2;

            % Create EmissionwavelengthEditFieldLabel
            app.EmissionwavelengthEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.EmissionwavelengthEditFieldLabel.HorizontalAlignment = 'right';
            app.EmissionwavelengthEditFieldLabel.Position = [18 539 118 22];
            app.EmissionwavelengthEditFieldLabel.Text = 'Emission wavelength';

            % Create EmissionwavelengthEditField
            app.EmissionwavelengthEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.EmissionwavelengthEditField.Limits = [100 Inf];
            app.EmissionwavelengthEditField.ValueDisplayFormat = '%11.4gnm';
            app.EmissionwavelengthEditField.ValueChangedFcn = createCallbackFcn(app, @EmissionwavelengthEditFieldValueChanged, true);
            app.EmissionwavelengthEditField.Position = [145 539 58 22];
            app.EmissionwavelengthEditField.Value = 100;

            % Create BeaddiameterEditFieldLabel
            app.BeaddiameterEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.BeaddiameterEditFieldLabel.HorizontalAlignment = 'right';
            app.BeaddiameterEditFieldLabel.Position = [53 626 83 22];
            app.BeaddiameterEditFieldLabel.Text = 'Bead diameter';

            % Create BeaddiameterEditField
            app.BeaddiameterEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.BeaddiameterEditField.ValueDisplayFormat = '%11.4gµm';
            app.BeaddiameterEditField.ValueChangedFcn = createCallbackFcn(app, @BeaddiameterEditFieldValueChanged, true);
            app.BeaddiameterEditField.Position = [145 626 58 22];
            app.BeaddiameterEditField.Value = 0.17;

            % Create RIsamplelayerEditFieldLabel
            app.RIsamplelayerEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.RIsamplelayerEditFieldLabel.HorizontalAlignment = 'right';
            app.RIsamplelayerEditFieldLabel.Position = [48 568 88 22];
            app.RIsamplelayerEditFieldLabel.Text = 'RI sample layer';

            % Create RIsamplelayerEditField
            app.RIsamplelayerEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.RIsamplelayerEditField.ValueChangedFcn = createCallbackFcn(app, @RIsamplelayerEditFieldValueChanged, true);
            app.RIsamplelayerEditField.Position = [145 568 58 22];
            app.RIsamplelayerEditField.Value = 1.33;

            % Create MagnificationEditFieldLabel
            app.MagnificationEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.MagnificationEditFieldLabel.HorizontalAlignment = 'right';
            app.MagnificationEditFieldLabel.Position = [282 626 76 22];
            app.MagnificationEditFieldLabel.Text = 'Magnification';

            % Create MagnificationEditField
            app.MagnificationEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.MagnificationEditField.Limits = [1 Inf];
            app.MagnificationEditField.ValueDisplayFormat = '%11.4gx';
            app.MagnificationEditField.ValueChangedFcn = createCallbackFcn(app, @MagnificationEditFieldValueChanged, true);
            app.MagnificationEditField.Position = [367 626 57 22];
            app.MagnificationEditField.Value = 100;

            % Create NAEditFieldLabel
            app.NAEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.NAEditFieldLabel.HorizontalAlignment = 'right';
            app.NAEditFieldLabel.Position = [333 597 25 22];
            app.NAEditFieldLabel.Text = 'NA';

            % Create NAEditField
            app.NAEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.NAEditField.Limits = [0.1 Inf];
            app.NAEditField.ValueChangedFcn = createCallbackFcn(app, @NAEditFieldValueChanged, true);
            app.NAEditField.Position = [367 597 57 22];
            app.NAEditField.Value = 1.2;

            % Create RIimmersionEditFieldLabel
            app.RIimmersionEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.RIimmersionEditFieldLabel.HorizontalAlignment = 'right';
            app.RIimmersionEditFieldLabel.Position = [282 568 76 22];
            app.RIimmersionEditFieldLabel.Text = 'RI immersion';

            % Create RIimmersionEditField
            app.RIimmersionEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.RIimmersionEditField.Limits = [0.01 Inf];
            app.RIimmersionEditField.ValueChangedFcn = createCallbackFcn(app, @RIimmersionEditFieldValueChanged, true);
            app.RIimmersionEditField.Position = [367 568 57 22];
            app.RIimmersionEditField.Value = 1.5;

            % Create PixelsizephysicalLabel
            app.PixelsizephysicalLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.PixelsizephysicalLabel.HorizontalAlignment = 'right';
            app.PixelsizephysicalLabel.Position = [248 539 110 22];
            app.PixelsizephysicalLabel.Text = 'Pixel size (physical)';

            % Create PixelsizephysicalEditField
            app.PixelsizephysicalEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.PixelsizephysicalEditField.Limits = [0.01 Inf];
            app.PixelsizephysicalEditField.ValueDisplayFormat = '%11.4gµm';
            app.PixelsizephysicalEditField.ValueChangedFcn = createCallbackFcn(app, @PixelsizephysicalEditFieldValueChanged, true);
            app.PixelsizephysicalEditField.Position = [367 539 57 22];
            app.PixelsizephysicalEditField.Value = 6.5;

            % Create FitresultsLabel
            app.FitresultsLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.FitresultsLabel.HorizontalAlignment = 'center';
            app.FitresultsLabel.FontSize = 14;
            app.FitresultsLabel.FontWeight = 'bold';
            app.FitresultsLabel.Visible = 'off';
            app.FitresultsLabel.Position = [686 473 72 22];
            app.FitresultsLabel.Text = 'Fit results';

            % Create LoadparametersButton
            app.LoadparametersButton = uibutton(app.PSFcharacterizationUIFigure, 'push');
            app.LoadparametersButton.ButtonPushedFcn = createCallbackFcn(app, @LoadparametersButtonValueChanged, true);
            app.LoadparametersButton.Position = [549 625 110 23];
            app.LoadparametersButton.Text = 'Load parameters';

            % Create parameterFile
            app.parameterFile = uieditfield(app.PSFcharacterizationUIFigure, 'text');
            app.parameterFile.Editable = 'off';
            app.parameterFile.Position = [667 626 227 22];

            % Create NumberZernikesEditFieldLabel
            app.NumberZernikesEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.NumberZernikesEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberZernikesEditFieldLabel.Position = [561 539 97 22];
            app.NumberZernikesEditFieldLabel.Text = 'Number Zernikes';

            % Create NumberZernikesEditField
            app.NumberZernikesEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.NumberZernikesEditField.Limits = [1 Inf];
            app.NumberZernikesEditField.RoundFractionalValues = 'on';
            app.NumberZernikesEditField.ValueChangedFcn = createCallbackFcn(app, @NumberZernikesEditFieldValueChanged, true);
            app.NumberZernikesEditField.Position = [667 539 32 22];
            app.NumberZernikesEditField.Value = 21;

            % Create PloterrorButton
            app.PloterrorButton = uibutton(app.PSFcharacterizationUIFigure, 'state');
            app.PloterrorButton.ValueChangedFcn = createCallbackFcn(app, @PloterrorButtonValueChanged, true);
            app.PloterrorButton.Enable = 'off';
            app.PloterrorButton.Visible = 'off';
            app.PloterrorButton.Text = 'Plot error';
            app.PloterrorButton.Position = [820 568 75 23];

            % Create SetcolorlimitspersliceCheckBox
            app.SetcolorlimitspersliceCheckBox = uicheckbox(app.PSFcharacterizationUIFigure);
            app.SetcolorlimitspersliceCheckBox.ValueChangedFcn = createCallbackFcn(app, @SetcolorlimitspersliceCheckBoxValueChanged, true);
            app.SetcolorlimitspersliceCheckBox.Text = 'Set color limits per slice';
            app.SetcolorlimitspersliceCheckBox.Position = [747 674 148 22];

            % Show the figure after all components are created
            app.PSFcharacterizationUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = aberration_measurement_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.PSFcharacterizationUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.PSFcharacterizationUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PSFcharacterizationUIFigure)
        end
    end
end