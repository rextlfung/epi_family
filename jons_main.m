%% create and plot SMS EPI sequence files
addpath ~/code/HarmonizedMRI/SMS-EPI/sequence/Pulseq/   % getsmspulse.m

writeSMSEPI2;

smsepi2ge;

system('tar xf epi.tar');
%toppe.plotseq(sysGE, 'timeRange', [0 0.06]);  % time range in sec
toppe.plotseq(sysGE);


%% create 3D GRE scan
write3DGRE;
gre3d2ge;

%% load data and interpolate onto Cartesian grid along x
fn = 'P,rex,epi.7';

loaddata;   % dat = [nx ny nSlices nCoils]
[nx ny nSlices nCoils] = size(dat);


%% ghost correction
d = squeeze(dat(:,:,3,:));  % ghost calibration data (phase encode off)
x = fftshift(ifft(fftshift(d), [], 1));  % getoephase expects image space
[a, th] = hmriutils.epi.getoephase(x, false);
dat_gc = hmriutils.epi.epiphasecorrect(dat, a);

%% recon and display
for sl = 1:size(dat_gc,3)
    figure;
    [~, I] = toppe.utils.ift3(squeeze(dat_gc(:,:,sl,:)), 'type', '2d');
    im(I)
end
