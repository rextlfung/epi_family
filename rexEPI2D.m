%% Simple low-performance 2D EPI sequence. Does not consider TE/TR (For educational purposes)
% Rex, Oct 11, 2023
clear; close all;

%% Parameters
fov = 21.6e-2; % 21.6 cm x 21.6 cm slices
Nx = 90; Ny = Nx; % 90 pixels x 90 pixels slices
sliceThickness = 2.4e-3; % 2.4 mm thick slices
Nslices = 5; % 5 slices in total
Nframes = 4; % 4 frames in total
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
          
%% 90 degree slice-selective excitation pulse
[rf, gz] = mr.makeSincPulse(tipAngle*pi/180,...
    'system',lims,...
    'Slicethickness',sliceThickness,...
    'Duration',rfDur,... 
    'timeBwProduct',rfTBW,...
    'apodization',0.5);

%% Other gradients and ADC events
readoutTime = Nx*dwellTime; % time that ADC is sampling
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',kWidth/readoutTime,'FlatTime',flatTime); % Readout gradient
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2); % ADC event

%% Pre-phasing gradients in x and y, and Re-phasing gradient in z
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2); % Rephase the dephasing from slice selection
dur = gzReph.riseTime + gzReph.flatTime + gzReph.fallTime;

gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2,'Duration',dur); % removed -deltak/2 to aligh the echo between the samples
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',dur); % Move to edge of k-space (-k_y,max)

%% Add Gy blips to move to the next line in k-space AFAP
grt = lims.gradRasterTime;
dur = ceil(2*sqrt(deltak/lims.maxSlew)/grt)*grt; % Shortest possible blip duration
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur); % Move to next line in k-sapce

%% String together sequence
seq = mr.Sequence(); % Create sequence object

for nframe = 1:Nframes
    
    for slice=1:Nslices 
        segmentID = slice; % Needed for seq2ge
        
        rf.freqOffset=gz.amplitude*sliceThickness*(slice-1-(Nslices-1)/2); % Shift RF frequency to shift slice location
        seq.addBlock(rf,gz,mr.makeLabel('SET', 'LIN', segmentID)); % Slice-selective excitation
        seq.addBlock(gxPre,gyPre,gzReph); % Prephase in x and y, rephase in z
        
        for i=1:Ny
            seq.addBlock(gx,adc);           % Read one line of k-space
            seq.addBlock(gy);               % Phase blip
            
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end
    end
    
    % Calculate additional delay needed to have desired TR
    if nframe == 1
        delayNeededForTR = TR - sum(seq.blockDurations);
        assert(delayNeededForTR >= 0);
    end
    seq.addBlock(mr.makeDelay(delayNeededForTR)); % Repetition time
end

seq.plot(); % Plot sequence waveforms

%% Calculate and plot k-space trajectory
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();

time_axis=(1:(size(ktraj,2)))*lims.gradRasterTime;
figure; plot(ktraj(1,:),ktraj(2,:),'b');
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.');

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
seq.write('rexEPI,2D.seq');

%% Convert to GE
sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
    'maxSlew', 20, ...                        % G/cm/ms
    'maxRF', 0.15, ...                  % Gauss. Must be >= peak RF in sequence.
    'maxView', Ny, ...               % Determines slice/view index in data file
    'maxSlice', Nslices*Nframes,...  % Saves the data processing hassle later
    'adcDeadTime', 20, ...           % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ...          % RF/gradient delay (us)
    'psd_grd_wait', 156);            % ADC/gradient delay (us)

seq2ge('rexEPI,2D.seq',sysGE,'rexEPI,2D.tar');

system('tar -xvf rexEPI,2D.tar');

figure; toppe.plotseq(sysGE,'timeRange',[0 sum(seq.blockDurations)/Nframes]); % Plot 1 frame