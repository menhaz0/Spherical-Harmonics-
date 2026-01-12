clear;

%grid
theta = linspace(0, pi, 100);
phi = linspace(0, 2*pi, 100);
[PHI, THETA] = meshgrid(phi, theta);

% --many surfaces in spherical coordinates -- 
%g = sin(THETA).^2 .* cos(2*PHI);
%g = 1 ./ (sin(THETA) + 0.1);
%g = cos(THETA)*sin(PHI);
%g = cos(THETA); - paraboloid
%g = sin(THETA + PHI); 
%g = double(THETA < pi/2) - double(THETA > pi/2); %Gibbs phenomenon for
%spherical coordinates? 

g = cos(THETA) ./ (sin(THETA).^2 + 0.1); %change g here

% -- calculing the coefficients using inner product and the complete set -- 

L = 20; %if you want less terms, change this variable. L = 20 could take a long time. 
A_lm = zeros(L + 1, 2*L + 1);
for l = 0:L

    for m = -l:l
        Y = spherical_harmonic(l, m, THETA, PHI);
        Norm = (4*pi / (2*l + 1)) * (factorial(l + abs(m)) / factorial(l - abs(m)));
        A_lm(l+1, m + L + 1) = (1/Norm) * trapz(phi, trapz(theta, g .* conj(Y) .* sin(THETA), 1));
    end
end

% -- plotting and summing -- 

g_lm = zeros(size(THETA));
fig = figure('Color', 'White', 'Position', [50, 50, 1200, 400]);
for l = 0:L
    for m = -l:l
        cf = A_lm(l+1, m + L + 1);
        y = spherical_harmonic(l, m, THETA, PHI);
        
        g_lm = g_lm + (cf * y);
        
        clf;  

        subplot(1,3,1);
        my_plot(THETA, PHI, g);
        title('$g(\theta, \phi)$', "Interpreter", "latex");
       
        subplot(1,3,2);
        my_plot(THETA, PHI, real(cf * y));
        title(sprintf('$A_{%d,%d} Y_{%d}^{%d}$ ($A_{%d,%d} \\approx %.4f$)',l,m,l,m,l,m,real(cf)),'Interpreter','latex');

        subplot(1,3,3);
        my_plot(THETA, PHI, real(g_lm)); %takes the real part of g_lm. 
        title(sprintf('Approximation: $g_{%d,%d}(\\theta, \\phi)$', l, m), 'Interpreter', 'latex');
        
        drawnow;
        pause(1); %changing how fast you want to iterate over the plots. 
    end
end

function my_plot(THETA, PHI, g)
    % spherical to rectangular coordinates 
    r = abs(g);
    X = r .* sin(THETA) .* cos(PHI);
    Y = r .* sin(THETA) .* sin(PHI);
    Z = r .* cos(THETA);

    surf(X, Y, Z, g);
    axis equal;
    grid on;

    axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]); %you can change this if you want to see g over a larger domain
end

% -- functions that return spherical harmonics evaluated at theta and phi--  
function Y = spherical_harmonic(l, m, theta, phi) %takes in inputs l, m, theta, and phi
    % l,m Legendre function evalulated at cos(theta)
    Plm = legendreFunction(l, m, cos(theta));
    Y = Plm.*exp(1i * m * phi);
end

% -- function that returns Legendre function evaluated at s-- 
function P = legendreFunction(l, m, s) %takes in inputs l, m, and s
    abs_m = abs(m);
    p0 = [1 0 -1]; %this vector represents s^2 + 0s - 1

    p_to_l = 1; % initializes to a constant polynomial

    for i = 1:l
        p_to_l = conv(p_to_l, p0); %this is just polynomial multiplication, but faster
    end

    p_der = p_to_l;
    %taking (l+m) derivatives
    for i = 1:(l + abs_m)
        p_der = polyder(p_der);
    end

    P = ((-1)^abs_m) / (2^l * factorial(l)).*(1 - s.^2).^(abs_m/2).* polyval(p_der, s);
end

