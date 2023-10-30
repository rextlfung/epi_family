%% Simple low-performance 3D EPI sequence. Does not consider TE/TR (For educational purposes)
% Rex, Oct 20, 2023
% v2 additions: water-selective excitation, 
clear; close all;

%% Parameters (ABCD)
Nframes = 4; % 4 repeated acquisitions
Nx = 90; Ny = Nx; Nz = 60; % 90 x 90 x 60 voxels
del_x = 2.4e-3; del_y = del_x; del_z = del_x; % 2.4 mm resolution isotropic
fov_x = Nx*del_x; fov_y = Ny*del_y; fov_z = Nz*del_z; % 216 mm x 216 mm x 144 mm volume
TR = 800e-3; % repetition time (s)
TE = 30e-3; % echo time (s)

% k-space params
del_kx =1/fov_x; del_ky = 1/fov_y; del_kz = 1/fov_z; % width of each voxel in k-sapce (1/m)
width_kx = Nx*del_kx; width_ky = Ny*del_ky; width_kz = Nz*del_kz; % toal width in k-space (1/m)

% other params
dwellTime = 16e-6; % ADC sample interval (s)
tipAngle = 52; % degrees
rfDur = 0.15e-3; % RF duration. Makes the amplitude just under 0.125 G for a 52 degree tip angle

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
          
%% volume excitation pulse, water selective
[rf, rfDelay] = mr.makeBlockPulse(tipAngle*pi/180,...
    'system',lims,...
    'Duration', rfDur);

waterFatFreqOffset = 440; % Hz (Derived from 3.5 ppm at 3T)
rfGap = 1/waterFatFreqOffset; % s (delay between first and second hard pulse)
rf.signal = rf.signal ./ 2; % half the RF amplitude

%% Other gradients and ADC events
readoutTime = Nx*dwellTime; % time that ADC is sampling
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',width_kx/readoutTime,'FlatTime',flatTime); % Readout gradient
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2); % ADC event

%% Pre-readout gradients in x and y, and Re-phasing gradient in z
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

%% Make all gradients GE compatible
crt = 20e-6; % Common raster time (multiple) of 10 us (Siemens) and 4 us (GE)
gxPre = trap4ge(gxPre,crt,lims);
gyPre = trap4ge(gyPre,crt,lims);
gzPre = trap4ge(gzPre,crt,lims);
gyBlip = trap4ge(gyBlip,crt,lims);

%% String together sequence
seq = mr.Sequence(); % Create sequence object

% Begin acquisition
for nframe = 1:Nframes
    % Loop through each plane
    for nz = 1:Nz 
        segmentID = nz; % Needed for seq2ge
        
        % Water-selective volume excitation
        seq.addBlock(rf,mr.makeDelay(rfGap),mr.makeLabel('SET', 'LIN', segmentID));
        seq.addBlock(rf);
        
        gzPreScaled = mr.scaleGrad(gzPre,nz - 1 - Nz/2); % Scale to proper z prephasing
        seq.addBlock(gxPre,gyPre,gzPreScaled); % PE in x, y, z
        
        for ny=1:Ny % EPI sample each plane
            seq.addBlock(gx,adc); % Read one line of k-space
            seq.addBlock(gyBlip); % Phase blip to next line
            
            gx.amplitude = -gx.amplitude; % Reverse polarity of read gradient
        end
    end
    
%     % Calculate additional delay needed to have desired TR
%     if nframe == 1
%         delayNeededForTR = TR - sum(seq.blockDurations);
%         assert(delayNeededForTR >= 0);
%     end
%     seq.addBlock(mr.makeDelay(delayNeededForTR)); % Repetition time
end

seq.plot('timeRange',[0, sum(seq.blockDurations)/Nframes]); % Plot sequence waveforms

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
seq.write('rexEPI,3D,v2.seq');

%% Convert to GE
sysGE = toppe.systemspecs('maxGrad', 5, ... % G/cm
    'maxSlew', 20, ... % G/cm/ms
    'maxRF', 0.15, ... % Gauss. Must be >= peak RF in sequence
    'maxView', Ny, ... % Determines slice/view index in data file
    'maxSlice', Nz*Nframes,... % Saves the data processing hassle later
    'adcDeadTime', 20, ... % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ... % RF/gradient delay (us)
    'psd_grd_wait', 156); % ADC/gradient delay (us)

seq2ge('rexEPI,3D,v2.seq',sysGE,'rexEPI,3D,v2.tar');

system('tar -xvf rexEPI,3D,v2.tar');

figure; toppe.plotseq(sysGE,'timeRange',[0, sum(seq.blockDurations)/Nframes]); % Plot 1 frame