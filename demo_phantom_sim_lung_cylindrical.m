%% Create a phantom
omega = 0.1;
osc_amp1 = 0.05;
osc_amp2 = 0.05;
osc_amp3 = 0.05;
osc_amp4 = 0.1;
osc_amp5 = 0.1;
osc_amp6 = 0.1;
tr = 0.015; % in s
phase = 0 * pi/180; % in radians
base_value = 1;
imSz = [128, 128, 128];
reconSz = [128, 128, 128];
reconstruct_flag = false;
mat_file = 'data_evolve_nonuniform_6cylinders_trajrand_nFrames_10000.mat';
%% Generate the trajectory
npts = 128;
dwell_time = 15;
ramp_time = 100;
plat_time = 2500;
decay_time = 60;
oversampling = 3;
del_x = 0;
del_y = 0;
del_z = 0;
flag = 3; % haltonized spiral
nFrames = 10000;
m_lNumberOfFrames = 1;

radialDistance_x = repelem(linspace(0, 1, npts), 3, 1)';
radialDistance_y = repelem(linspace(0, 1, npts), 3, 1)';
radialDistance_z = repelem(linspace(0, 1, npts), 3, 1)';

[x_end, y_end, z_end] = GX_f_gen_traj(nFrames, m_lNumberOfFrames, flag);
traj_scale_factor = 1;
x = traj_scale_factor * 0.5 * imSz(1) * radialDistance_x(:, 1) * x_end;
y = traj_scale_factor * 0.5 * imSz(2) * radialDistance_y(:, 2) * y_end;
z = traj_scale_factor * 0.5 * imSz(3) * radialDistance_z(:, 3) * z_end;

traj = reshape(z/imSz(1), [npts * nFrames, 1]);
traj(:, 2) = reshape(x/imSz(2), [npts * nFrames, 1]);
traj(:, 3) = reshape(y/imSz(3), [npts * nFrames, 1]);

%% Sample k-space
kspIm = zeros(npts, nFrames);
for i = 1:nFrames
    intensity1 = 2 * base_value + osc_amp1 * square(i*omega + phase);
    intensity2 = base_value + osc_amp2 * square(i*omega + phase);
    intensity3 = base_value + osc_amp3 * square(i*omega + phase);
    intensity4 = base_value + osc_amp4 * square(i*omega + phase);
    intensity5 = base_value + osc_amp5 * square(i*omega + phase);
    intensity6 = base_value + osc_amp6 * square(i*omega + phase);
    st = mri_objects('cyl3', [-.15, 0, -0.2, 0.1, 0.2, intensity1], ...
        'cyl3', [-0.15, 0, 0, 0.1, 0.2, intensity2], ...
        'cyl3', [-0.15, 0, 0.2, 0.1, 0.2, intensity3], ...
        'cyl3', [0.15, 0, .2, 0.1, 0.2, intensity4], ...
        'cyl3', [0.15, 0, 0, 0.1, 0.2, intensity5], ...
        'cyl3', [0.15, 0, -.2, 0.1, 0.2, intensity6]);
    kspIm_i = st.kspace(x(:, i), y(:, i), z(:, i));
    kspIm(:, i) = squeeze(kspIm_i);
end
data = reshape(kspIm, [npts * nFrames, 1]);

%% Save mat file
[X, Y, Z] = meshgrid(linspace(-0.5, .5, imSz(1)), ...
    linspace(-0.5, .5, imSz(2)), ...
    linspace(-0.5, .5, imSz(3)));
image = st.image(X, Y, Z);
mask = imrotate3(double(image > osc_amp1 * 0.1), 90, [1, 0, 0]);
st = mri_objects('cyl3', [-.15, 0, -0.2, 0.1, 0.2, osc_amp1], ...
    'cyl3', [-0.15, 0, 0, 0.1, 0.2, osc_amp2], ...
    'cyl3', [-0.15, 0, 0.2, 0.1, 0.2, osc_amp3], ...
    'cyl3', [0.15, 0, .2, 0.1, 0.2, osc_amp4], ...
    'cyl3', [0.15, 0, 0, 0.1, 0.2, osc_amp5], ...
    'cyl3', [0.15, 0, -.2, 0.1, 0.2, osc_amp6]);
image = st.image(X, Y, Z);
mask_evolve = imrotate3(image, 90, [1, 0, 0]);
save(mat_file, ...
    'traj', ...
    'data', ...
    'mask', ...
    'mask_evolve');

%% Reconstruction
if reconstruct_flag
    kernel.sharpness = 0.32;
    kernel.extent = 9 * kernel.sharpness;
    overgrid_factor = 3;
    nDcfIter = 15;
    deapodizeImage = false; %true();
    cropOvergriddedImage = true();
    verbose = true();

    %  Choose kernel, proximity object, and then create system model
    kernelObj = Recon.SysModel.Kernel.Gaussian(kernel.sharpness, kernel.extent, verbose);
    proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);
    clear kernelObj;

    systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
        reconSz, proxObj, verbose);

    % Choose density compensation function (DCF)
    dcfObj = Recon.DCF.Iterative(systemObj, nDcfIter, verbose);

    % Choose Reconstruction Model
    reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
    clear modelObj;
    clear dcfObj;
    reconObj.crop = cropOvergriddedImage;
    reconObj.deapodize = deapodizeImage;

    % Reconstruct image using trajectories in pixel units
    recon = reconObj.reconstruct(data, traj);

    %% Visualization
    figure(1);
    imslice(image);
end

%% More Visualization
close all;
figure;
frame = 1;
plot3([-64, 64], [0, 0], [0, 0], 'k', 'LineWidth', 2); hold on;
plot3([0, 0], [-64, 64], [0, 0], 'k','LineWidth', 2);
plot3([0, 0], [0, 0], [-64, 64], 'k', 'LineWidth', 2);
while frame <= 25
        if square(frame*omega) > 0
            plot3([0, x(end, frame)], [0, y(end, frame)], [0, z(end, frame)], 'k-*');
        else
            plot3([0, x(end, frame)], [0, y(end, frame)], [0, z(end, frame)], 'k-*');
        end
    frame = frame + 1;
    % pause(0.1);
    hold on;
end
pbaspect([1 1 1]);
figure;
frame = 1;
while frame <= nFrames
    if frame <= nFrames / 2
        plot3([0, x(end, frame)], [0, y(end, frame)], [0, z(end, frame)], 'k.');
    else
        plot3([0, x(end, frame)], [0, y(end, frame)], [0, z(end, frame)], 'r.');
    end
    frame = frame + 1;
    % pause(0.2);
    hold on;
end

%% Calculate the sum of distance between ends of the radial trajectories
% interleaved
x1 = x(end, 1:2:end);
x2 = x(end, 2:2:end);
y1 = y(end, 1:2:end);
y2 = y(end, 2:2:end);
z1 = z(end, 1:2:end);
z2 = z(end, 2:2:end);
distance_matrix = zeros(nFrames/2, nFrames/2);
for i = 1:nFrames/2
    for j = 1:nFrames/2
        distance_matrix(i, j) = sqrt((x1(i)-x2(j))^2+(y1(i)-y2(j))^2+(z1(i)-z2(j))^2);
    end
end

disp(mean(min(distance_matrix, [], 2)))

% repeated
x1 = x(end, 1:nFrames/2);
x2 = x(end, nFrames/2:end);
y1 = y(end, 1:nFrames/2);
y2 = y(end, nFrames/2:end);
z1 = z(end, 1:nFrames/2);
z2 = z(end, nFrames/2:end);
distance_matrix = zeros(nFrames/2, nFrames/2);
for i = 1:nFrames/2
    for j = 1:nFrames/2
        distance_matrix(i, j) = sqrt((x1(i)-x2(j))^2+(y1(i)-y2(j))^2+(z1(i)-z2(j))^2);
    end
end

disp(mean(min(distance_matrix, [], 2)))