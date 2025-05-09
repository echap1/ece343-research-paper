% --- Grid and constants ---
Nx = 1000; Ny = 1000;
dx = 1; dy = 1;
c0 = 299792458;
eps0 = 8.854e-12; 
mu0 = 4*pi*1e-7;

dt = 1 / (c0 * sqrt((1/dx^2)+(1/dy^2))) * 0.99;  % CFL limit
Nt = 1000;

% --- Fields ---
Ez = zeros(Nx, Ny);
Hx = zeros(Nx, Ny-1);
Hy = zeros(Nx-1, Ny);

% --- Constants ---
CHx = dt / (mu0 * dy);
CHy = dt / (mu0 * dx);
CEz = dt / eps0;

% --- Source ---
src_x = round(Nx/2);
src_y = round(Ny/2);

Ez_max = 0.005;        % Approx. peak of Gaussian (exp(0) = 1)
H_max = Ez_max / sqrt(mu0 / eps0);  % From wave impedance: |E| / |H| = η = sqrt(mu0 / eps0)


% --- Object Setup: Dielectric square in upper left ---
eps_r = ones(Nx, Ny);  % Relative permittivity
block_size = 200;
eps_r(100:100+block_size, 100:100+block_size) = 4.0;
eps_r(100:100+block_size, end-100-block_size:end-100) = 2.0;

% Update CEz to reflect material variation
CEz_matrix = dt ./ (eps0 * eps_r);

% --- PEC circle in bottom right ---
[X, Y] = meshgrid(1:Ny, 1:Nx);
circle_center = [800, 800];
circle_radius = 50;
pec_mask = (X - circle_center(2)).^2 + (Y - circle_center(1)).^2 < circle_radius^2;


% --- Visualization of Simulation Space ---
figure('Name', 'Simulation Geometry', 'Position', [100 100 800 800]);

% Visualize relative permittivity
imagesc(eps_r');
axis equal tight;
% colormap('hot');
colorbar;
title('Simulation Geometry: Dielectric and PEC Objects');
hold on;

% Overlay PEC circle
theta = linspace(0, 2*pi, 200);
x_circ = circle_center(2) + circle_radius * cos(theta);
y_circ = circle_center(1) + circle_radius * sin(theta);
plot(x_circ, y_circ, 'k-', 'LineWidth', 2);

legend('PEC boundary');


figure('Position', [100 100 800 800]);

gif_filename = 'objects.gif';


% --- Time loop ---
for n = 1:Nt
    % Update magnetic fields
    Hx(:,1:end) = Hx(:,1:end) - CHx * diff(Ez,1,2);  % dEz/dy
    Hy(1:end,:) = Hy(1:end,:) + CHy * diff(Ez,1,1);  % dEz/dx

    % Calculate curl(H)
    curl_H = zeros(Nx, Ny);
    curl_H(2:end-1,2:end-1) = diff(Hy(:,2:end-1),1,1) - diff(Hx(2:end-1,:),1,2);

    % Update electric field
    Ez = Ez + CEz_matrix .* curl_H;
    Ez(pec_mask) = 0;

    % % --- First-order Mur ABC on boundaries ---
    % if n > 1
    %     % Top and bottom
    %     Ez(1,:)   = Ez_old(2,:)   + ((c0*dt - dx)/(c0*dt + dx)) * (Ez(2,:)   - Ez_old(1,:));
    %     Ez(end,:) = Ez_old(end-1,:) + ((c0*dt - dx)/(c0*dt + dx)) * (Ez(end-1,:) - Ez_old(end,:));
    % 
    %     % Left and right
    %     Ez(:,1)   = Ez_old(:,2)   + ((c0*dt - dy)/(c0*dt + dy)) * (Ez(:,2)   - Ez_old(:,1));
    %     Ez(:,end) = Ez_old(:,end-1) + ((c0*dt - dy)/(c0*dt + dy)) * (Ez(:,end-1) - Ez_old(:,end));
    % end

    % Store current Ez for Mur ABC use in next step
    Ez_old = Ez;

    % Source
    Ez(src_x, src_y) = Ez(src_x, src_y) + exp(-((n-30)/10)^2);

    % Plot every 5 steps
    if mod(n,5) == 0
        % Resize H fields to match Ez dimensions
        Hx_plot = padarray(Hx, [0 1], 'post');  % Nx x Ny
        Hy_plot = padarray(Hy, [1 0], 'post');  % Nx x Ny
    
        % Compute |H| magnitude
        H_mag = sqrt(Hx_plot.^2 + Hy_plot.^2);
    
        % Fixed normalization using physical expected max
        Ez_clipped = max(min(abs(Ez), Ez_max), 0);  % Clamp for display
        H_mag_clipped = min(H_mag, H_max);           % Clamp to avoid oversaturation
        
        R = Ez_clipped / Ez_max;    % Map [0, Ez_max] to [0,1]
        B = H_mag_clipped / H_max;                   % Map [0, H_max] to [0,1]
        G = zeros(size(R));

        % Combine into RGB image
        RGB = cat(3, R', G', B');  % transpose to match display
    
        image(RGB); 
        axis equal tight off;
        title(['RGB field at timestep ', num2str(n)]);
        drawnow;
    
        % Save frame to GIF
        frame = getframe(gcf);
        im = frame2im(frame);
        im = imresize(im, 0.25);  % Scale to 25% of original size (adjust as needed)
        [imind, cm] = rgb2ind(im, 256);
        if n == 5  % First frame (corresponds to mod(n,5)==0 and n==5)
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
        end
    end

end
