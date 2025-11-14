% recon and display 2D EPI multislice images
% Assumes readout module is 'readout.mod'

% Experimental parameters
fn = 'P,epi2d,Nz=10,16Sep2023.7';
Nz = 10;
Nx = 64;
% Name readout module 'readout.mod'

% load data and interpolate onto Cartesian grid along x
loaddata;   % dat = [nx ny nSlices nCoils]
[nx ny nSlices nCoils] = size(dat);

clear I
for sl = 1:Nz
    % ghost correction for this slice
    d = squeeze(dat(:,:,sl,:));  % ghost calibration data (phase encode off)
    x = fftshift(ifft(fftshift(d), [], 1));  % getoephase expects image space
    verbose = false;
    [a, th] = hmriutils.epi.getoephase(x, verbose);
    datc = hmriutils.epi.epiphasecorrect(squeeze(dat(:,:,sl+Nz,:)), a);

    % recon and display
    [~, I(:,:,sl)] = toppe.utils.ift3(datc, 'type', '2d');
end

I = flipdim(flipdim(I,1),2);  % to match orientation on scanner host
figure; im(I);
