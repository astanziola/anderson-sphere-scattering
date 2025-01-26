function P = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, order, options)
    % COMPUTEANDERSONSPHERESOLUTION Computes the acoustic scattering from a fluid sphere
    %   This function calculates the analytical solution for acoustic scattering
    %   from a fluid sphere in a fluid medium using Anderson's formulation.
    %   The solution assumes that the sphere is centered at the origin and that
    %   the incident wave is a plane wave propagating in the negative z-direction.
    %   The field can only be evaluated outside the sphere.
    %
    % Arguments:
    %   p: Position matrix [3 x N] where each column is a point [x; y; z] (m)
    %   c0: Speed of sound in the surrounding medium (m/s)
    %   rho0: Density of the surrounding medium (kg/m³)
    %   c1: Speed of sound in the sphere (m/s)
    %   rho1: Density of the sphere (kg/m³)
    %   R: Radius of the sphere (m)
    %   omega: Angular frequency of the incident wave (rad/s)
    %   order: Maximum order of spherical harmonics to include
    
    arguments
        p (:,:) {mustBeNumeric}
        c0 {mustBeNumeric, mustBePositive}
        rho0 {mustBeNumeric, mustBePositive}
        c1 {mustBeNumeric, mustBePositive}
        rho1 {mustBeNumeric, mustBePositive}
        R {mustBeNumeric, mustBePositive}
        omega {mustBeNumeric, mustBePositive}
        order {mustBeNumeric, mustBeInteger, mustBeNonnegative}
        options.kind {mustBeMember(options.kind, {'plane-wave', 'point-source'})} = 'plane-wave'
        options.D {mustBeNumeric, mustBePositive} = 1.0  % Point source distance on the z-axis
    end
    
    % Ensure p is 3xN
    if size(p, 1) ~= 3
        error('Position matrix p must have 3 rows [x; y; z]');
    end
    
    % Convert Cartesian to spherical coordinates (vectorized)
    r = sqrt(sum(p.^2, 1));  % Compute radius for each point
    theta = acos(p(3,:)./r);  % Compute theta for each point
    
    validatePosition(r, R);
    
    % Calculate wave properties (vectorized)
    k0 = omega/c0;                     % Wavenumber in medium
    k1 = omega/c1;                     % Wavenumber in sphere
    kr = k0*r;                         % Dimensionless radius
    mu = cos(theta);                   % Directivity variable
    
    % Compute scattered field using modal expansion
    use_point_source = strcmp(options.kind, 'point-source');
    P = computeModalExpansion(order, k0, k1, R, rho0, c0, rho1, c1, kr, mu, use_point_source, options.D);
    
    % Add incident wave
    P = P + computeIncidentWave(k0, r, mu, use_point_source, options.D);
end

function P0 = computeIncidentWave(k0, r, mu, use_point_source, D)
    % Compute the incident wave using a plane wave or point source (vectorized)
    if ~use_point_source
        P0 = exp(-1i .* k0 .* r .* mu);
    else
        % Monopole source located in (0, 0, D)
        r0 = sqrt(D.^2 + r.^2 - 2 .* D .* r .* mu);
        P0 = exp(1i .* k0 .* r0) ./ (1i .* k0 .* r0);
    end
end

function P = computeModalExpansion(order, k0, k1, R, rho0, c0, rho1, c1, kr, mu, use_point_source, D)
    % Compute the scattered field using modal expansion (vectorized)
    P = zeros(size(kr));  % Initialize output array
    
    for m = 0:order
        % Compute modal coefficient (scalar)
        Am = computeModalCoefficient(m, k0, k1, R, rho0, c0, rho1, c1, use_point_source, D);
        
        % Compute spherical wave functions (vectorized)
        h1_m = computeSphericalHankel(m, kr);
        Pm = evaluateLegendrePolynomial(m, mu);
        
        % Add modal contribution (vectorized)
        P = P + Am * Pm .* h1_m;
    end
end

function Am = computeModalCoefficient(m, k0, k1, r, rho0, c0, rho1, c1, use_point_source, D)
    % Modal coefficient computation remains the same as it's scalar
    if use_point_source
        Lm = computeSphericalHankel(m, k0.*D);
    else
        Lm = (-1i).^m;
    end

    Z0 = rho0 * c0;
    Z1 = rho1 * c1;
    
    k0r = k0*r;
    k1r = k1*r;
    
    j1_m_k1r = computeSphericalBesseljDerivative(m, k1r);
    jm_k0r = computeSphericalBesselj(m, k0r);
    j1_m_k0r = computeSphericalBesseljDerivative(m, k0r);
    jm_k1r = computeSphericalBesselj(m, k1r);
    hm_k0r = computeSphericalHankel(m, k0r);
    h1_m_k0r = computeSphericalHankelDerivative(m, k0r);
    
    numerator = j1_m_k1r * jm_k0r * Z0 - j1_m_k0r * jm_k1r * Z1;
    denominator = -j1_m_k1r * hm_k0r * Z0 + h1_m_k0r * jm_k1r * Z1;
    
    Am = (2*m + 1) * Lm * numerator / denominator;
end

function j1_m = computeSphericalBesseljDerivative(m, x)
    % Compute derivative of spherical Bessel function of the first kind
    % See: https://functions.wolfram.com/Bessel-TypeFunctions/BesselJ/20/ShowAll.html
    j1_m = 0.5 * (computeSphericalBesselj(m-1, x) - computeSphericalBesselj(m+1, x));
end

function y1_m = computeSphericalBesselyDerivative(m, x)
    % Compute derivative of spherical Bessel function of the second kind
    % See: https://functions.wolfram.com/Bessel-TypeFunctions/BesselY/20/ShowAll.html
    y1_m = 0.5 * (computeSphericalBessely(m-1, x) - computeSphericalBessely(m+1, x));
end

function h1_m = computeSphericalHankelDerivative(m, x)
    % Compute derivative of spherical Hankel function of the first kind
    h1_m = computeSphericalBesseljDerivative(m, x) + ...
           1i * computeSphericalBesselyDerivative(m, x);
end

function hm = computeSphericalHankel(m, x)
    % Compute spherical Hankel function of the first kind
    hm = computeSphericalBesselj(m, x) + 1i * computeSphericalBessely(m, x);
end

function jm = computeSphericalBesselj(m, x)
    % COMPUTESPHERICALBESSELJ Compute spherical Bessel function of the first kind
    %   See: https://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html
    
    validateattributes(x, {'numeric'}, {'finite', 'nonnan'}, ...
                      'computeSphericalBesselj', 'x');
    
    Jm = besselj(m + 0.5, x);
    jm = sqrt(pi./(2*x)) .* Jm;
end

function ym = computeSphericalBessely(m, x)
    % COMPUTESPHERICALBESSELJ Compute spherical Bessel function of the second kind
    %   Also known as the spherical Neumann function
    %   See: https://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
    
    validateattributes(x, {'numeric'}, {'finite', 'nonnan'}, ...
                      'computeSphericalBessely', 'x');
    
    Ym = bessely(m + 0.5, x);
    ym = sqrt(pi./(2*x)) .* Ym;
end

function validatePosition(r, R)
    % Validate that all evaluation points are outside the sphere (vectorized)
    if any(r <= R)
        error('AndersonSolution:InvalidPosition', ...
              'All evaluation points must be outside the sphere (r > R)');
    end
end

% Custom validation functions
function mustBeInteger(x)
    if any(mod(x, 1) ~= 0)
        error('Value must be an integer');
    end
end

function mustBePositive(x)
    if any(x <= 0)
        error('Value must be positive');
    end
end