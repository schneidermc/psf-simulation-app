classdef aberration_measurement_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PSFcharacterizationUIFigure   matlab.ui.Figure
        NumberZernikesEditField       matlab.ui.control.NumericEditField
        NumberZernikesEditFieldLabel  matlab.ui.control.Label
        parameterFile                 matlab.ui.control.EditField
        LoadparametersButton          matlab.ui.control.StateButton
        FitresultsLabel               matlab.ui.control.Label
        PixelsizephysicalEditField    matlab.ui.control.NumericEditField
        PixelsizephysicalLabel        matlab.ui.control.Label
        RIimmersionEditField          matlab.ui.control.NumericEditField
        RIimmersionEditFieldLabel     matlab.ui.control.Label
        NAEditField                   matlab.ui.control.NumericEditField
        NAEditFieldLabel              matlab.ui.control.Label
        MagnificationEditField        matlab.ui.control.NumericEditField
        MagnificationEditFieldLabel   matlab.ui.control.Label
        RIfluidEditField              matlab.ui.control.NumericEditField
        RIfluidEditFieldLabel         matlab.ui.control.Label
        BeaddiameterEditField         matlab.ui.control.NumericEditField
        BeaddiameterEditFieldLabel    matlab.ui.control.Label
        EmissionwavelengthEditField   matlab.ui.control.NumericEditField
        EmissionwavelengthEditFieldLabel  matlab.ui.control.Label
        zincrementEditField           matlab.ui.control.NumericEditField
        zincrementEditFieldLabel      matlab.ui.control.Label
        IterationsEditField           matlab.ui.control.NumericEditField
        MaxiterationsLabel            matlab.ui.control.Label
        CalibrationsampleLabel        matlab.ui.control.Label
        MicroscopeparametersLabel     matlab.ui.control.Label
        text_stackfile                matlab.ui.control.EditField
        ErrorLabel                    matlab.ui.control.Label
        Lamp                          matlab.ui.control.Lamp
        LoadzstackButton              matlab.ui.control.StateButton
        FitButton                     matlab.ui.control.Button
        zSliderSimulation             matlab.ui.control.Slider
        CalculatemodelButton          matlab.ui.control.Button
        zSliderMeasurement            matlab.ui.control.Slider
        FitZernikeaberrationsLabel    matlab.ui.control.Label
        SaveparametersButton          matlab.ui.control.Button
        SaveaberrationsButton         matlab.ui.control.Button
        SavetransmissionButton        matlab.ui.control.Button
        UIAxes_Transmission           matlab.ui.control.UIAxes
        UIAxes_Zernike                matlab.ui.control.UIAxes
        UIAxes_PhaseAberration        matlab.ui.control.UIAxes
        UIAxesSimulationPsfXY         matlab.ui.control.UIAxes
        UIAxesSimulationPsfXZ         matlab.ui.control.UIAxes
        UIAxesMeasurementPsfXZ        matlab.ui.control.UIAxes
        UIAxesMeasurementPsfXY        matlab.ui.control.UIAxes
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
        PSF_xz % max. intensity projection
        PSF_xy
        PSF_simu_xz
        PSF_simu_xy
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
        function display_projection(app, axes, data_xz, Nx, Ny, Nz, ux, dz)
            x_axis = linspace(0, (Nx-1)*ux*1e6, Nx);
            y_axis = linspace(0, (Ny-1)*ux*1e6, Ny);
            z_axis = linspace(0, (Nz-1)*dz*1e6, Nz);
            % Show xz projection
            h = axes; % select figure
            imagesc(h, x_axis, z_axis, data_xz);
            set(h,'DataAspectRatio',[dz, ux, ux])
            axis(h, "tight");
            colormap(h,app.CallingApp.ColormapDropDown.Value)
        end
               
        function display_single_image(app, axes, Nx, Ny, ux, data)
            h = axes; % select figure
            x_axis = linspace(0, (Nx-1)*ux*1e6, Nx);
            y_axis = linspace(0, (Ny-1)*ux*1e6, Ny);
            imagesc(h,x_axis, y_axis, data);
            set(h,'DataAspectRatio',[ux, ux, ux]);
            axis(h,"tight");
            colormap(h,app.CallingApp.ColormapDropDown.Value)
        end
        
        function drawZLineMeasurement(app,value)
            persistent lineHandleMeasurement
            if isempty(lineHandleMeasurement) || ~isvalid(lineHandleMeasurement)
                lineHandleMeasurement = yline(app.UIAxesMeasurementPsfXZ, value,'w-');
            end
            % Update line position
            lineHandleMeasurement.Value = value; 
        end

        function drawZLineSimulation(app,value)
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
            
            app.RIfluidEditField.Value = app.CallingApp.RefractiveIndexSpecimenSpinner.Value;
            app.RI_fluid = app.RIfluidEditField.Value;
            
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

        % Value changed function: RIfluidEditField
        function RIfluidEditFieldValueChanged(app, event)
            app.RI_fluid = app.RIfluidEditField.Value;
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

        % Value changed function: LoadparametersButton
        function LoadparametersButtonValueChanged(app, event)
            [file, path] = uigetfile('*.mat');

            if isempty(who('-file', fullfile(path,file), "parameters"))
                errordlg('No parameters found in selected file!','Load error');
            else
                load(fullfile(path,file),"parameters");
                app.BeaddiameterEditField.Value = parameters.beadDiameter;
                app.zincrementEditField.Value = parameters.zIncrement;
                app.RIfluidEditField.Value = parameters.RiFluid;
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

        % Value changed function: LoadzstackButton
        function LoadzstackButtonValueChanged(app, event)
            % Loading stack of bead images
            [filename, pathname] = uigetfile('*.tif');
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
            app.PSF_xz = rot90(squeeze(sum(app.stack,2))); %standard projection
            
            % Make items visible
            app.UIAxesMeasurementPsfXZ.Visible = "on";
            app.UIAxesMeasurementPsfXY.Visible = "on";
            %app.objectiveclosesttobeadLabel.Visible = "on";
            %app.objectivefarthestfrombeadLabel.Visible = "on";

            % Display the bead images
            display_projection(app, app.UIAxesMeasurementPsfXZ, app.PSF_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, app.stack(:,:,round(app.zSliderMeasurement.Value)));
            
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
            app.PSF_simu_xz = rot90(squeeze(sum(app.stack_simu,2))); % standard projection
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
                show_error(app, app.stack, app.stack_simu);
                app.FitButton.Enable = "on";
            else
                val = app.Nz;
            end
            flippedValue = app.Nz + 1 - val;
            app.zSliderSimulation.Value = val;

            display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
            display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, app.stack_simu(:,:,val));
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
            outputFunction = @(x,optimValues,state) updateProgress(x,optimValues,state,d,resnorm_initial);
            options = optimoptions(@lsqnonlin,'Algorithm',algorithm ,'MaxFunctionEvaluations', app.iter_max,'OutputFcn',outputFunction);

            function stop = updateProgress(x,optimValues,state,d,resnorm_initial)
                % Update the Value property of the progress bar to reflect the progress
                % based on the decrease in the residual norm
                resnorm_current = optimValues.resnorm;
                if resnorm_current < resnorm_initial
                    progress = (resnorm_initial - resnorm_current) / resnorm_initial;
                else
                    progress = 0.01;
                end
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
            
            % Zernikes
            stem(app.UIAxes_Zernike, app.modes(3:end), app.Z_phase(3:end)); % omit 2nd and 3rd (tip/tilt) as they are set to 0
            grid(app.UIAxes_Zernike, "on");
            
            % Check retrieved bead image
            [~, app.stack_simu,~]=fun_errorcalc(app.Z_aberr, app.stack/sum(app.stack(:)), alpha, apo_no, kernel, F_kernel, defocus, z_vec, app.E_BFP, app.dia_pupil, Zernike_stack, F_bead, app.slice_norm);
            
            % Display simulated bead images with aberrations taken into account
            app.PSF_simu_xz = rot90(squeeze(sum(app.stack_simu,2))); %standard projection
            display_projection(app, app.UIAxesSimulationPsfXZ, app.PSF_simu_xz, app.Nx, app.Ny, app.Nz, app.ux, app.dz);
            display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, app.stack_simu(:,:,1));
            app.zSliderSimulation.Value = app.Nz;
            drawZLineMeasurement(app, 0);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, app.stack(:,:,1));
            app.zSliderMeasurement.Value = app.Nz;
            drawZLineSimulation(app, 0);
            
            show_error(app, app.stack, app.stack_simu);
            app.SaveaberrationsButton.Enable = "on";
            app.SavetransmissionButton.Enable = "on";
        end

        % Value changed function: zSliderSimulation
        function zSliderSimulationValueChanged(app, event)
            app.zSliderSimulation.Value = round(event.Value);
            app.zSliderMeasurement.Value = app.zSliderSimulation.Value;
            flippedValue = app.Nz + 1 - app.zSliderSimulation.Value;
            data = app.stack_simu(:,:,flippedValue);
            display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, data)
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, app.stack(:,:,flippedValue))
            app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
            app.drawZLineSimulation((flippedValue-1)*app.zincrementEditField.Value)
        end

        % Value changing function: zSliderSimulation
        function zSliderSimulationValueChanging(app, event)
            app.zSliderSimulation.Value = round(event.Value);
            app.zSliderMeasurement.Value = app.zSliderSimulation.Value;
            flippedValue = app.Nz + 1 - app.zSliderSimulation.Value;
            data = app.stack_simu(:,:,flippedValue);
            display_single_image(app, app.UIAxesSimulationPsfXY, app.Nx, app.Ny, app.ux, data)
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, app.stack(:,:,flippedValue))
            app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
            app.drawZLineSimulation((flippedValue-1)*app.zincrementEditField.Value)
        end

        % Value changed function: zSliderMeasurement
        function zSliderMeasurementValueChanged(app, event)
            app.zSliderMeasurement.Value = round(event.Value);
            flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
            data = app.stack(:,:,flippedValue);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, data)
            app.drawZLineMeasurement((flippedValue-1)*app.zincrementEditField.Value)
        end

        % Value changing function: zSliderMeasurement
        function zSliderMeasurementValueChanging(app, event)
            app.zSliderMeasurement.Value = round(event.Value);
            flippedValue = app.Nz + 1 - app.zSliderMeasurement.Value;
            data = app.stack(:,:,flippedValue);
            display_single_image(app, app.UIAxesMeasurementPsfXY, app.Nx, app.Ny, app.ux, data)
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
            transmission = app.transmissionMask;            
            [file,path] = uiputfile('*.mat','Save results','transmission.mat');
            save(fullfile(path,file),"transmission")
        end

        % Value changed function: NumberZernikesEditField
        function NumberZernikesEditFieldValueChanged(app, event)
            app.NumberZernikesEditField.Value = event.Value;
            app.modes = 2:app.NumberZernikesEditField.Value;
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
            app.UIAxes_Zernike.Box = 'on';
            app.UIAxes_Zernike.YGrid = 'on';
            app.UIAxes_Zernike.Visible = 'off';
            app.UIAxes_Zernike.Position = [549 96 345 185];

            % Create UIAxes_Transmission
            app.UIAxes_Transmission = uiaxes(app.PSFcharacterizationUIFigure);
            title(app.UIAxes_Transmission, 'Transmission')
            app.UIAxes_Transmission.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_Transmission.XTick = [];
            app.UIAxes_Transmission.YTick = [];
            app.UIAxes_Transmission.ZTick = [];
            app.UIAxes_Transmission.Box = 'on';
            app.UIAxes_Transmission.Visible = 'off';
            colormap(app.UIAxes_Transmission, 'gray')
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
            app.LoadzstackButton = uibutton(app.PSFcharacterizationUIFigure, 'state');
            app.LoadzstackButton.ValueChangedFcn = createCallbackFcn(app, @LoadzstackButtonValueChanged, true);
            app.LoadzstackButton.VerticalAlignment = 'top';
            app.LoadzstackButton.Text = 'Load z-stack';
            app.LoadzstackButton.Position = [548 597 110 23];

            % Create Lamp
            app.Lamp = uilamp(app.PSFcharacterizationUIFigure);
            app.Lamp.Visible = 'off';
            app.Lamp.Position = [874 541 20 20];
            app.Lamp.Color = [0.651 0.651 0.651];

            % Create ErrorLabel
            app.ErrorLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.ErrorLabel.Visible = 'off';
            app.ErrorLabel.Position = [667 568 154 22];
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

            % Create RIfluidEditFieldLabel
            app.RIfluidEditFieldLabel = uilabel(app.PSFcharacterizationUIFigure);
            app.RIfluidEditFieldLabel.HorizontalAlignment = 'right';
            app.RIfluidEditFieldLabel.Position = [93 568 43 22];
            app.RIfluidEditFieldLabel.Text = 'RI fluid';

            % Create RIfluidEditField
            app.RIfluidEditField = uieditfield(app.PSFcharacterizationUIFigure, 'numeric');
            app.RIfluidEditField.ValueChangedFcn = createCallbackFcn(app, @RIfluidEditFieldValueChanged, true);
            app.RIfluidEditField.Position = [145 568 58 22];
            app.RIfluidEditField.Value = 1.33;

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
            app.LoadparametersButton = uibutton(app.PSFcharacterizationUIFigure, 'state');
            app.LoadparametersButton.ValueChangedFcn = createCallbackFcn(app, @LoadparametersButtonValueChanged, true);
            app.LoadparametersButton.VerticalAlignment = 'top';
            app.LoadparametersButton.Text = 'Load parameters';
            app.LoadparametersButton.Position = [548 626 110 23];

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