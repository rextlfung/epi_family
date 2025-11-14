% SMS-EPI sequence in Pulseq
%
% This script creates the file 'smsepi.seq', that can be executed directly
% on Siemens MRI scanners using the Pulseq interpreter.
% The .seq file can also be converted to a .tar file that can be executed on GE
% scanners, see main.m.
%
% The experimental parameters below are chosen such that the sequence 
% can be executed identically (to us precision) on Siemens and GE systems.
% For more information about preparing a Pulseq file for execution on GE scanners,
% see the 'Pulseq on GE' manual.

scanType = 'SMS';   % 'SMS' or '3D'

warning('OFF', 'mr:restoreShape');

%% Paths
caipiPythonPath = '~/code/HarmonizedMRI/3DEPI/caipi/';

%% Define experimental parameters
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 120, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 20e-6, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 10e-6, ...
              'blockDurationRaster', 10e-6, ...
              'B0', 3.0);

% ABCD:
% voxelSize = [2.4 2.4 2.4]*1e-3; 
% mb = 6;
% Nx = 90; Ny = Nx; Nz = mb*10;
% TE = 30e-3;
% alpha = 52;
voxelSize = [3.5 3.5 3.5]*1e-3;     % m
mb = 1;                             % multiband/SMS factor
Nx = 64; Ny = Nx;                   % Matrix size
Nz = 10*mb; 
TE = 30e-3;                         % echo time (s)
alpha = 52;                         % flip angle (degrees)

fov = voxelSize .* [Nx Ny Nz];      % FOV (m)
slThick = fov(3)/Nz

nFrames = 2;                       % number of temporal frames (image volumes)
nDummyFrames = 0;                   % dummy frames to reach steady state

dwell = 4e-6;                       % ADC sample time (s). For GE, must be multiple of 2us.

% CAIPI sampling parameters
Ry = 1;                   % ky undersampling factor
Rz = mb;                  % kz undersampling factor
CaipiShiftZ = 2;
pf_ky = 1.0;              % partial Fourier factor along ky
etl = ceil(pf_ky*Ny/Ry);  % echo train length

nCyclesSpoil = 2;    % number of spoiler cycles, along x and z

rfTB  = 6;          % RF pulse time-bandwidth product
rfDur = 8e-3;       % RF pulse duration (s)

fatChemShift = 3.5*1e-6;                        % 3.5 ppm
fatOffresFreq = sys.gamma*sys.B0*fatChemShift;  % Hz


%% Create fat sat pulse 

% TODO


%% Create excitation pulse
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
    'maxRF', 0.25);

if strcmp(scanType, 'SMS')
    sliceSep = fov(3)/mb   % center-to-center separation between SMS slices (m)
    % freq = frequency offset (Hz) corresponding to a sliceSep
    [rf, gzRF, freq] = getsmspulse(alpha, slThick, rfTB, rfDur, ...
        mb, sliceSep, sysGE, sys, ...
        'doSim', true, ...    % Plot simulated SMS slice profile
        'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
        'ftype', 'ls');       % filter design. 'ls' = least squares

    freq = freq/sliceSep*voxelSize(3); % frequency offset for z shift of voxelSize(3)
end


%% Get CAIPI sampling pattern (for one shot/echo train)

% create caipi.mat, and load it
pyFile = [caipiPythonPath 'skippedcaipi_sampling.py'];
pyCmd = sprintf('python3 %s %d %d %d %d %d %d', ...
    pyFile, Ny, Nz, Ry, Rz, CaipiShiftZ, 1);
a = input('Press 1 to run Python script from Matlab, 2 to run offline:  ');
if a == 1
    system(pyCmd);
else
    fprintf('Open a terminal and run the following python command:\n\t%s\n', pyCmd);
    input('Then press Enter to continue');
end

load caipi

% kz and ky indeces (multiples of deltak)
kyInds = double(indices((end-etl+1):end, 2));
kzInds = double(indices((end-etl+1):end, 1));

% ky/kz encoding blip amplitude along echo train (multiples of deltak)
kyStep = diff(kyInds);
kzStep = diff(kzInds);

kyStepMax = max(abs(kyStep));
kzStepMax = max(abs(kzStep));

%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

deltak = 1./fov;

% Start with the blips
gyBlip = mr.makeTrapezoid('y', sys, 'Area', kyStepMax*deltak(2)); 
gzBlip = mr.makeTrapezoid('z', sys, 'Area', kzStepMax*deltak(3)); 

% Area and duration of the biggest blip
if gyBlip.area > gzBlip.area
    maxBlipArea = gyBlip.area;
    blipDuration = mr.calcDuration(gyBlip);
else
    maxBlipArea = gzBlip.area;
    blipDuration = mr.calcDuration(gzBlip);
end

% Readout trapezoid
systmp = sys;
systmp.maxGrad = deltak(1)/dwell;  % to ensure >= Nyquist sampling
gro = mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea);

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
if mod(round(Tread*1e6)*1e-6, dwell)
    Tread = Tread - mod(Tread, dwell) + dwell;
end
adc = mr.makeAdc(round(Tread/dwell), sys, ...
    'Duration', Tread, ...
    'Delay', blipDuration/2);

% Split blips at block boundary
[gyBlipUp, gyBlipDown] = mr.splitGradientAt(gyBlip, mr.calcDuration(gyBlip)/2);
gyBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gyBlipUp);
gyBlipDown.delay = 0;

[gzBlipUp, gzBlipDown] = mr.splitGradientAt(gzBlip, mr.calcDuration(gzBlip)/2);
gzBlipUp.delay = mr.calcDuration(gro) - mr.calcDuration(gzBlipDown);
gzBlipDown.delay = 0;

% prephasers and spoilers
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gro.area/2);
Tpre = mr.calcDuration(gxPre);
gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', (kyInds(1)-Ny/2)*deltak(2), ... 
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', -Nz/2*deltak(3), ... 
    'Duration', Tpre);
%gxSpoil = mr.makeTrapezoid('x', sys, ...
%    'Area', Nx*deltak(1)*nCyclesSpoil);
gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);


%% Calculate delay to achieve desired TE
kyIndAtTE = find(kyInds-Ny/2 == min(abs(kyInds-Ny/2)));
minTE = mr.calcDuration(gzRF) - mr.calcDuration(rf)/2 - rf.delay + mr.calcDuration(gxPre) + ...
        (kyIndAtTE-0.5) * mr.calcDuration(gro);
TEdelay = floor((TE-minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;


%% Interleaved slice ordering for SMS
% if even number of shots, swap last two to reduce slice cross-talk
IZ = 1:Rz:Nz;
if strcmp(scanType, 'SMS')
    Nex = length(IZ);  % number of excitations to cover one image volume
    II = [1:2:Nex 2:2:Nex];
    if ~mod(Nex,2)
        l = length(II);
        II = II([1:(l-2) l l-1]);
    end

    % for testing:
    II = IZ;
else
    % 3D EPI
    II = 1:length(IZ);
end


%% Assemble sequence
seq = mr.Sequence(sys);           

nShots = Nz/mb;

% temporal frame loop
for ifr = 0:nFrames

    % dummy shots to reach steady state before receive gain calibration in Auto Prescan
    isDummyShot = ifr < 1;

    % If 2D EPI, acquire EPI reference scan (turn off gy blips)
    isRefShot = ifr == 1 & mb == 1;

    % Segment ID (see 'Pulseq on GE' manual)
    segmentID = 3;
    if isDummyShot
        segmentID = 1;
    end
    if isRefShot
        segmentID = 2;
    end

    blipsOn = ~isDummyShot & ~isRefShot;

    % slice (group) offset loop
    for ii = II

        if strcmp(scanType, 'SMS')
            rf.freqOffset = round((ii-0.5-length(IZ)/2)*freq);  % frequency offset (Hz) for slice shift
            a_gz = 0;  % prephaser amplitude
        else
            rf.freqOffset = 0;
            a_gz = (1-isRefShot)*((IZ(ii)-1+0.5)-Nz/2)/(Nz/2);
        end

        % RF excitation
        % TODO: RF spoiling
        seq.addBlock(rf, gzRF, mr.makeLabel('SET', 'LIN', segmentID));

        % TE delay
        if TE > minTE
            seq.addBlock(mr.makeDelay(TEdelay));
        end

        % Readout
        %seq.addBlock(gxPre, gyPre, mr.scaleGrad(gzPre, 1 - 2/Nz*(shot-1)));
        seq.addBlock(gxPre, mr.scaleGrad(gyPre, blipsOn), mr.scaleGrad(gzPre, strcmp(scanType,'3D')*isRefShot));

        if isDummyShot
            seq.addBlock(gro, ...
                         mr.scaleGrad(gyBlipUp, blipsOn*kyStep(1)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipUp, blipsOn*kzStep(1)/max(kzStepMax,1)));
        else
            seq.addBlock(gro, adc, ...
                         mr.scaleGrad(gyBlipUp, blipsOn*kyStep(1)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipUp, blipsOn*kzStep(1)/max(kzStepMax,1)));
        end

        for ie = 2:(etl-1)
            % Example: 'gybdu' = Gy blip down up
            gybd = mr.scaleGrad(gyBlipDown, blipsOn*kyStep(ie-1)/max(kyStepMax,1));
            gybu = mr.scaleGrad(gyBlipUp,   blipsOn*kyStep(ie)/max(kyStepMax,1));
            gybdu = mr.addGradients({gybd, gybu}, sys);
            gzbd = mr.scaleGrad(gzBlipDown, blipsOn*kzStep(ie-1)/max(kzStepMax,1));
            gzbu = mr.scaleGrad(gzBlipUp,   blipsOn*kzStep(ie)/max(kzStepMax,1));
            gzbdu = mr.addGradients({gzbd, gzbu}, sys);
            if isDummyShot
                seq.addBlock(mr.scaleGrad(gro, (-1)^(ie-1)), gybdu, gzbdu);
            else
                seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(ie-1)), gybdu, gzbdu);
            end
        end

        if isDummyShot
            seq.addBlock(mr.scaleGrad(gro, (-1)^(ie)), ...
                         mr.scaleGrad(gyBlipDown, blipsOn*kyStep(ie)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipDown, blipsOn*kzStep(ie)/max(kzStepMax,1)));
        else
            seq.addBlock(adc, ...
                         mr.scaleGrad(gro, (-1)^(ie)), ...
                         mr.scaleGrad(gyBlipDown, blipsOn*kyStep(ie)/max(kyStepMax,1)), ...
                         mr.scaleGrad(gzBlipDown, blipsOn*kzStep(ie)/max(kzStepMax,1)));
        end

        % spoil. Comment out to make the k-space plot neater
        seq.addBlock(gzSpoil);

        % Long TR for testing
        %seq.addBlock(mr.makeDelay(1.0));
    end
 end

%% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Output for execution and plot
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'smsepi');
seq.write('smsepi.seq')       % Write to pulseq file

seq.plot('timeRange', [0 0.06]);

%% k-space trajectory calculation and plot
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

%% Optional slow step, but useful for testing during development,
%% e.g., for the real TE, TR or for staying within slewrate limits
% rep = seq.testReport;
% fprintf([rep{:}]);
