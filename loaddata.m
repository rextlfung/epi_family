addpath ~/github/HarmonizedMRI/utils/

fn = '/mnt/storage/rexfung/epi/P,rexEPI3D,231027.7';

% get raw data for 1 'slice' (for testing)
d = toppe.utils.loadpfile(fn);
d = flip(d, 1);      % [nfid ncoils nslices nechoes nviews]
d = permute(d, [1 5 3 2 4]);   % [nfid Ny nslices nCoils nechoes]
d = squeeze(d);  % [nfid nCoils nslices nviews]

size(d)

%% Discard calibration data
d = d(:,:,2:end,:);
[Nx Ny Nz Ncoils] = size(d)

%% Split z and time dimensions
Nframes = 2;
d = reshape(d,[Nx Ny Nz/Nframes Nframes Ncoils]);

%% Flip every other line in ky
d(:,2:2:end,:,:,:) = flip(d(:,2:2:end,:,:,:), 1);

%% Check max k-space value for clipping
d_vec = d(:);
fprintf('max real value: %d\n',max(real(d_vec)));
fprintf('fraction of data points near clipping: %d\n',...
    sum(real(d_vec) > 32700)/length(d_vec));

%% Look at k-space
figure;
for coil = 1:Ncoils
    subplot(1,2,1);
    im(abs(d(:,:,:,1,coil)).^0.3);
    subplot(1,2,2);
    im(abs(d(:,:,:,2,coil)).^0.3);
    pause;
end

%% Look at image
img = toppe.utils.ift3(squeeze(d(:,:,:,1,:)));
img_ss = sum(conj(img).*img,4).^0.5;
figure; im(img_ss)