%% Simple low-performance 3D EPI sequence. Does not consider TE/TR (For educational purposes)
% Rex, Oct 17, 2023
% Last edited Oct 27, 2023
clear; close all;

%% Parameters (ABCD)
Nframes = 2; % 4 repeated acquisitions
Nx = 90; Ny = Nx; Nz = 60; % 90 x 90 x 60 voxels
del_x = 2.4e-3; del_y = del_x; del_z = del_x; % 2.4 mm resolution isotropic
fov_x = Nx*del_x; fov_y = Ny*del_y; fov_z = Nz*del_z; % 216 mm x 216 mm x 144 mm volume
TR = 4.8; % repetition time (s)
TE = 30e-3; % echo time (s)

% k-space params
del_kx =1/fov_x; del_ky = 1/fov_y; del_kz = 1/fov_z; % width of each voxel in k-sapce (1/m)
width_kx = Nx*del_kx; width_ky = Ny*del_ky; width_kz = Nz*del_kz; % toal width in k-space (1/m)

% other params
Ndummy = 2; % Number of dummy scans (volumes) to drive magnetization to steady state
Ncal = 90; % Number of calibration lines to tune receiver gain
dwellTime = 4e-6; % ADC sample interval (s)
tipAngle = 52; % degrees
rfDur = 0.3e-3; % RF duration. Makes the amplitude just under 0.125 G for a 52 degree tip angle

lims = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
    'maxSlew', 110, 'slewUnit', 'T/m/s', ...
    'rfDeadTime', 100e-6, ...
    'rfRingdownTime', 60e-6, ...
    'adcDeadTime', 20e-6, ...
    'adcRasterTime', 2e-6, ...
    'gradRasterTime', 10e-6, ...
    'blockDurationRaster', 10e-6, ...
    'B0', 3.0);
          
% impose stronger slew constraint for spoilers since they are more likely to cause PNS
lims_spoiler = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
    'maxSlew', 60, 'slewUnit', 'T/m/s', ...
    'rfDeadTime', lims.rfDeadTime, ...
    'rfRingdownTime', lims.rfRingdownTime, ...
    'adcDeadTime', lims.adcDeadTime, ...
    'gradRasterTime', lims.gradRasterTime, ...
    'blockDurationRaster', lims.blockDurationRaster, ...
    'B0', lims.B0);
          
%% volume excitation pulse
[rf, rfDelay] = mr.makeBlockPulse(tipAngle*pi/180,...
    'system',lims,...
    'Duration', rfDur);

%% Other gradients and ADC events
readoutTime = Nx*dwellTime; % time that ADC is sampling
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',width_kx/readoutTime,'FlatTime',flatTime); % Readout gradient
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2); % ADC event

%% Pre-readout gradients in x and y, and Re-phasing gradient in z
% gzReph = mr.makeTrapezoid('z',lims,'Area',-gzEx.area/2); % Rephase the dephasing from slab selection
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2); % PE in x i.e. move to -k_x_max
gxDur = gxPre.riseTime + gxPre.flatTime + gxPre.fallTime;

% PE in y, i.e. move -k_y_max
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*del_ky,'Duration',gxDur);

% PE in z
gzPre = mr.makeTrapezoid('z',lims,'Area',-Nz/2*del_kz,'Duration',gxDur);
gzPre = mr.scaleGrad(gzPre,-2/Nz); % Gradient needed to move by 1 delta_k_z
%% Add Gy blips to move to the next line in k-space AFAP
grt = lims.gradRasterTime;
dur = ceil(2*sqrt(del_ky/lims.maxSlew)/grt)*grt; % Shortest possible blip duration
gyBlip = mr.makeTrapezoid('y',lims,'Area',del_ky,'Duration',dur); % Blip to next line in k-space

%% Spoiler gradients (move to -4*kx_max, 4*ky_max)
gxSpoil = mr.makeTrapezoid('x',lims_spoiler,'Area',3*-Nx/2*del_kx);
gySpoil = mr.makeTrapezoid('y',lims_spoiler,'Area',3*Ny/2*del_ky);
gzSpoil = mr.makeTrapezoid('z',lims_spoiler,'Area',3*Ny/2*del_kz);

%% Make all gradients GE compatible
crt = 20e-6; % Common raster time (multiple) of 10 us (Siemens) and 4 us (GE)
% gzReph = trap4ge(gzReph,crt,lims);
gxPre = trap4ge(gxPre,crt,lims);
gyPre = trap4ge(gyPre,crt,lims);
gzPre = trap4ge(gzPre,crt,lims);
gyBlip = trap4ge(gyBlip,crt,lims);
gxSpoil = trap4ge(gxSpoil,crt,lims_spoiler);
gySpoil = trap4ge(gySpoil,crt,lims_spoiler);

%% String together sequence
seq = mr.Sequence(); % Create sequence object

% % Dummy excitations to drive magnetization to steady state
% segmentID = 1; % needed for seq2ge
% for nframe = 1:Ndummy
%    seq.addBlock(rf,mr.makeLabel('SET', 'LIN', segmentID),mr.makeDelay(90e-3));
% end
% 
% % Acquire once with no PE for receiver gain tuning
% segmentID = 2;
% seq.addBlock(rf,mr.makeLabel('SET', 'LIN', segmentID));
% for line = 1:Ncal
%     seq.addBlock(adc); % Read a 90 x 90 plane of the center of k-space
% end
% %seq.addBlock(mr.makeDelay(87e-3)); % Let magnetization recover until actual data acquisition

% Begin acquisition
segmentID = 3;
for nframe = -Ndummy:Nframes
    % Loop through each plane
    for nz = 1:Nz 
        if nframe < 0 % Dummy excitations
            segmentID = 1;
        elseif nframe == 0 % Tune receiver gain
            segmentID = 2;
        else % Acquisition
            segmentID = 3; 
        end
        
        if nframe > 0
            % RF spoling
            % integer multiples of 117 degrees are usually indivisible by 360
            % quadratic phase cycling
            rf.phaseOffset = mod(117*(nz^2 + nz + 2)*pi/180, 2*pi);
            adc.phaseOffset = rf.phaseOffset;
        end
        
        seq.addBlock(rf,mr.makeLabel('SET', 'LIN', segmentID)); % Volume excitation
        
        gzPreScaled = mr.scaleGrad(gzPre,nz - 1 - Nz/2); % Scale to proper z prephasing
        if nframe > 0
            seq.addBlock(gxPre,gyPre,gzPreScaled); % PE in x, y, z
        else
            seq.addBlock(mr.scaleGrad(gxPre,0),...
                         mr.scaleGrad(gyPre,0),...
                         mr.scaleGrad(gzPreScaled,0));
        end
        
        if nframe > 0
            % EPI sample each kx-ky plane
            for ny=1:Ny
                seq.addBlock(gx,adc); % Read one line of k-space

                % Phase blip to next line, except for last line
                if ny < Ny 
                    seq.addBlock(gyBlip);
                end

                gx.amplitude = -gx.amplitude; % Reverse polarity of read gradient
            end
        elseif nframe == 0
            % Sample the center of k-space over and over again
            for ny=1:Ny
                seq.addBlock(mr.scaleGrad(gx,0),adc); % Read one line of k-space

                % Phase blip to next line, except for last line
                if ny < Ny 
                    seq.addBlock(mr.scaleGrad(gyBlip,0));
                end

                gx.amplitude = -gx.amplitude; % Reverse polarity of read gradient
            end
        else
            % Don't sample anything
            for ny=1:Ny
                seq.addBlock(mr.scaleGrad(gx,0)); % Read one line of k-space

                % Phase blip to next line, except for last line
                if ny < Ny 
                    seq.addBlock(mr.scaleGrad(gyBlip,0));
                end

                gx.amplitude = -gx.amplitude; % Reverse polarity of read gradient
            end
        end
        
        if nframe > 0
            % Spoiling in x and y and z
            seq.addBlock(gxSpoil,gySpoil,gzSpoil);
        else
            seq.addBlock(mr.scaleGrad(gxSpoil,0),...
                         mr.scaleGrad(gySpoil,0),...
                         mr.scaleGrad(gzSpoil,0));
        end
        
        % End z-loop if not actual acquisition
        if nframe <= 0
            break;
        end
    end
    
%     % Calculate additional delay needed to have desired TR
%     if nframe == 1
%         delayNeededForTR = TR - sum(seq.blockDurations);
%         assert(delayNeededForTR >= 0);
%     end
%     seq.addBlock(mr.makeDelay(delayNeededForTR)); % Repetition time
end

% seq.plot(); % Plot sequence waveforms

%% Calculate and plot k-space trajectory
if Ny <= 8
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();

    time_axis=(1:(size(ktraj,2)))*lims.gradRasterTime;
    figure; plot3(ktraj(1,:),ktraj(2,:),ktraj(3,:),'b');
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    hold; plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'r.');
    xlabel('k_x (cycles/m)'); ylabel('k_y (cycles/m)'); zlabel('k_z (cycles/m)');
end

%% Write to .seq file
seq.write('rexEPI,3D.seq');

%% Convert to GE
sysGE = toppe.systemspecs('maxGrad', 5, ... % G/cm
    'maxSlew', 20, ... % G/cm/ms
    'maxRF', 0.15, ... % Gauss. Must be >= peak RF in sequence
    'maxView', Ny, ... % Determines slice/view index in data file
    'maxSlice', 1 + Nz*Nframes,... % Saves the data processing hassle later
    'adcDeadTime', 20, ... % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ... % RF/gradient delay (us)
    'psd_grd_wait', 156); % ADC/gradient delay (us)

seq2ge('rexEPI,3D.seq',sysGE,'rexEPI,3D.tar');

system('tar -xvf rexEPI,3D.tar');

figure; toppe.plotseq(sysGE,'timeRange',[0, sum(seq.blockDurations)]); % Plot