% smsepi2ge.m
% Convert smsepi.seq to .tar file for execution on GE

sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...                          % G/cm/ms
    'maxRF', 0.15, ...                  % Gauss. Must be >= peak RF in sequence.
    'maxView', Ny, ...               % Determines slice/view index in data file
    'adcDeadTime', 20, ...           % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ...          % RF/gradient delay (us)
    'psd_grd_wait', 156);            % ADC/gradient delay (us)

seq2ge('smsepi.seq', sysGE, 'epi.tar');
