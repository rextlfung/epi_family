% recon and display 2D EPI multislice images
% Assumes readout module is 'readout.mod'

% Experimental parameters
fn = 'P,epi2d,mb=3,Nz=1,16Sep2023.7';
Nz = 1;
% Name readout module 'readout.mod'

% load data and interpolate onto Cartesian grid along x
loaddata;   % dat = [nx ny nSlices nCoils]
[nx ny nSlices nCoils] = size(dat)

% ghost correction
d = squeeze(dat(:,:,ceil(Nz/2),:));  % ghost calibration data (phase encode off)
x = fftshift(ifft(fftshift(d), [], 1));  % getoephase expects image space
verbose = true;
[a, th] = hmriutils.epi.getoephase(x, verbose);
datc = hmriutils.epi.epiphasecorrect(dat, a);

% recon and display
for sl = (Nz+1):(2*Nz)
    [~, I] = toppe.utils.ift3(squeeze(datc(:,:,sl,:)), 'type', '2d');
    I = flipdim(flipdim(I,1),2);  % to match orientation on scanner host
    figure; im(I); title(num2str(sl));
end
