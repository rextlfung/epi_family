%% Simple low-performance 3D EPI sequence. Does not consider TE/TR (For educational purposes)
% Rex, Oct 17, 2023
clear; close all;

%% Parameters
fov = 21.6e-2; % 21.6 cm x 21.6 cm x 21.6 cm volume
Nx = 16; Ny = Nx; Nz = Nx; % 64 x 64 x 64 voxels
Nframes = 4; % 4 volumes
deltak=1/fov; % width of each voxel in k-sapce (1/m)
kWidth = Nx*deltak; % toal width in k-space (1/m)
dwellTime = 4e-6; % ADC sample interval (s)
tipAngle = 90; % degrees
rfDur = 8e-3; % RF pulse duration (s)
rfTBW = 6; % RF time-bandwidth product
TR = 0.5; % Repetition time (s)

lims = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
    'maxSlew', 110, 'slewUnit', 'T/m/s', ...
    'rfDeadTime', 100e-6, ...
    'rfRingdownTime', 60e-6, ...
    'adcDeadTime', 20e-6, ...
    'adcRasterTime', 2e-6, ...
    'gradRasterTime', 10e-6, ...
    'blockDurationRaster', 10e-6, ...
    'B0', 3.0);
          
% impose stronger slew constraint for spoilers since they are more likely
% to cause PNS
lims_spoiler = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
    'maxSlew', 60, 'slewUnit', 'T/m/s', ...
    'rfDeadTime', lims.rfDeadTime, ...
    'rfRingdownTime', lims.rfRingdownTime, ...
    'adcDeadTime', lims.adcDeadTime, ...
    'gradRasterTime', lims.gradRasterTime, ...
    'blockDurationRaster', lims.blockDurationRaster, ...
    'B0', lims.B0);
          
%% 90 degree volume excitation pulse
[rf, gzEx] = mr.makeSincPulse(tipAngle*pi/180,...
    'system',lims,...
    'Slicethickness',fov,...
    'Duration',rfDur,... 
    'timeBwProduct',rfTBW,...
    'apodization',0.5);

%% Other gradients and ADC events
readoutTime = Nx*dwellTime; % time that ADC is sampling
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',kWidth/readoutTime,'FlatTime',flatTime); % Readout gradient
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2); % ADC event

%% Pre-phasing gradients in x and y, and Re-phasing gradient in z
gzReph = mr.makeTrapezoid('z',lims,'Area',-gzEx.area/2); % Rephase the dephasing from slab selection
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2); % Prephase in x i.e. move to -k_x_max

% Prephase in y an z, i.e. move to corner of 3D k-space
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak);
gzPre = mr.makeTrapezoid('z',lims,'Area',-Nz/2*deltak);

gzPre = mr.addGradients({gzReph,gzPre},'system',lims); % combine rephasing and prephasing gradients into one

%% Add Gy blips to move to the next line in k-space AFAP
grt = lims.gradRasterTime;
dur = ceil(2*sqrt(deltak/lims.maxSlew)/grt)*grt; % Shortest possible blip duration
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur); % Blip to next line in k-space
gzRo = mr.makeTrapezoid('z',lims,'Area',deltak,'Duration',dur); % Blip to next plane in k-space

%% String together sequence
seq = mr.Sequence(); % Create sequence object

segmentID = 1; % Needed for seq2ge

for nframe = 1:Nframes % Acquire a time series
    seq.addBlock(rf,gzEx,mr.makeLabel('SET', 'LIN', segmentID)); % Volume excitation
    seq.addBlock(gxPre,gyPre,gzPre); % Prephase in x, y, z. Rephase in z
    for nz = 1:Nz % Loop through each plane
        for ny=1:Ny % EPI sample each plane
            seq.addBlock(gx,adc);           % Read one line of k-space
            seq.addBlock(gy);               % Phase blip to next line
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end
        seq.addBlock(gzRo); % Phase blip to next plane
        gy.amplitude = -gy.amplitude; % Reverse polarity of Gy
    end
   seq.addBlock(mr.makeDelay(TR)); % Delay for T1 weighting, set to 0.5 s right now
end

seq.plot(); % Plot sequence waveforms

%% Calculate and plot k-space trajectory
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();

time_axis=(1:(size(ktraj,2)))*lims.gradRasterTime;
figure; plot3(ktraj(1,:),ktraj(2,:),ktraj(3,:),'b');
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'r.');
xlabel('k_x (cycles/m)'); ylabel('k_y (cycles/m)'); zlabel('k_z (cycles/m)');

%% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Write to .seq file
seq.write('rexEPI,3D,oneshot.seq');

%% Convert to GE
sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
    'maxSlew', 20, ...                        % G/cm/ms
    'maxRF', 0.15, ...                  % Gauss. Must be >= peak RF in sequence.
    'maxView', Ny, ...               % Determines slice/view index in data file
    'maxSlice', Nz*Nframes,...          % Saves the data processing hassle later
    'adcDeadTime', 20, ...           % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ...          % RF/gradient delay (us)
    'psd_grd_wait', 156);            % ADC/gradient delay (us)

seq2ge('rexEPI,3D,oneshot.seq',sysGE,'rexEPI,3D,oneshot.tar');

system('tar -xvf rexEPI,3D,oneshot.tar');

figure; toppe.plotseq(sysGE,'timeRange',[0, sum(seq.blockDurations)/Nframes - TR]); % Plot 1 frame