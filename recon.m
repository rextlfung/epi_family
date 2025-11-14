addpath ~/github/HarmonizedMRI/utils/

fn = '/mnt/storage/rexfung/epi/P,rexEPI3D,231027.7';

% get raw data for 1 'slice' (for testing)
rawdata = toppe.utils.loadpfile(fn);
rawdata = flip(rawdata, 1);      % [nfid ncoils nslices nechoes nviews]
rawdata = permute(rawdata, [1 5 3 2 4]);   % [nfid Ny nslices nCoils nechoes]
rawdata = squeeze(rawdata);  % [nfid nCoils nslices nviews]

fprintf('size of raw data: %d %d %d %d \n', size(rawdata))

%% Discard calibration data
d_cal = rawdata(:,:,1,:);
d = rawdata(:,:,2:end,:);
[Nx Ny NzNframes Ncoils] = size(d);

fprintf('size of data after discarding calibration data: %d %d %d %d \n', size(d))

%% Flip every other line in ky along kx, since it's zig-zag EPI
d(:,2:2:end,:,:,:) = flip(d(:,2:2:end,:,:,:), 1);

%% Check max k-space value for clipping
d_vec = d(:);
fprintf('max real value: %d\n',max(real(d_vec)));
fprintf('fraction of data points near clipping: %d\n',...
    sum(real(d_vec) > 32700)/length(d_vec));

%% Split z and time dimensions
Nz = 60;
Nframes = NzNframes/Nz;
d = reshape(d,[Nx Ny Nz Nframes Ncoils]);

fprintf('size of data after splitting z and time dimensions: %d %d %d %d %d \n', size(d))

%% Look at k-space
figure;
for coil = 16
    subplot(1,2,1);
    im(abs(d(:,:,:,1,coil)).^0.3);
    subplot(1,2,2);
    im(abs(d(:,:,:,2,coil)).^0.3);
end

%% Naive ghost correction: circshift the peak of each kx line to center
d_gc = d;
peakposes = zeros(1,Ncoils*Nframes*Nz*Ny);
for ncoil = 1:Ncoils
    disp(ncoil)
    for nframe = 1:Nframes
        for nz = 1:Nz
            for ny = 1:Ny
                mags = abs(d_gc(:,ny,nz,nframe,ncoil)); % Copy the current line of magnitudes along kx
                maxmag = max(mags); % Find the peak magnitude
                peakpos = round(median(find(mags == maxmag))); % Find peak position as median of all peaks (in case there's multiple)
                
                peakposes(ncoil*nframe*nz*ny) = peakpos;
                d_gc(:,ny,nz,nframe,ncoil) = circshift(d_gc(:,ny,nz,nframe,ncoil),Nx/2 - peakpos);
            end
        end
    end
end

% Verdict: FAILED. Enchanced ghosts even more.

%% Stats
datastats(peakposes')
histogram(peakposes)
%% Look at k-space
figure;
for coil = 1:Ncoils
    subplot(1,2,1);
    im(abs(d_gc(:,:,:,1,coil)).^0.3);
    subplot(1,2,2);
    im(abs(d_gc(:,:,:,2,coil)).^0.3);
    pause;
end

%% Look at image
img = toppe.utils.ift3(squeeze(d(:,:,:,2,:)));
img_ss = sum(conj(img).*img,4).^0.5;
figure; im(img_ss)