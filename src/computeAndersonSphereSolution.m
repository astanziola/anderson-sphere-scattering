function P = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, order, options)
    % COMPUTEANDERSONSPHERESOLUTION Computes the acoustic scattering from a fluid sphere
    %   This function calculates the analytical solution for acoustic scattering
    %   from a fluid sphere in a fluid medium using Anderson's formulation.
    %   The solution assumes that the sphere is centered at the origin and that
    %   the incident wave is a plane wave propagating in the negative z-direction.
    %   The field can only be evaluated outside the sphere.
    %
    % References:
    %   [1] Anderson V. C. (1950). "Sound scattering from a fluid sphere." The 
    %       Journal of the Acoustical Society of America, July 1950; 22 (4): 
    %       426–431
    %   [2] McNew  J., Lavarello R., O’Brien W. "Sound scattering from two 
    %       concentric fluid spheres." The Journal of the Acoustical Society 
    %       of America, 2007 Nov; 122(5):2968

    % Arguments:
    %   p: Position vector [x; y; z] where to evaluate the solution (m)
    %   c0: Speed of sound in the surrounding medium (m/s)
    %   rho0: Density of the surrounding medium (kg/m³)
    %   c1:  Speed of sound in the sphere (m/s)
    %   rho1: Density of the sphere (kg/m³)
    %   R: Radius of the sphere (m)
    %   omega: Angular frequency of the incident wave (rad/s)
    %   order: Maximum order of spherical harmonics to include
    %
    % Returns:
    %   P: Complex pressure amplitude at the specified position
    %
    % Notes:
    %   - The solution is valid for r > R, where r is the distance from the center
    %     of the sphere (origin) and R is the radius of the sphere.
    %
    % Example:
    %   p = [0.1; 0; 0.2];  % Position to evaluate
    %   c0 = 1.0;           % Water sound speed
    %   rho0 = 1.0;         % Water density
    %   c1 = 2.0;           % Sphere sound speed
    %   rho1 = 2.0;         % Sphere density
    %   R = 0.1;            % 5cm radius sphere
    %   f = 10;             % 10 Hz
    %   omega = 2*pi*f;     % Angular frequency
    %   order = 10;         % Order of the expansion
    %   P = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, order);
    %
    % ABOUT:
    %     author       - Antonio Stanziola
    %     date         - 11th January 2025
    %     last update  - 11th January 2025
    
    arguments
        p (3,1) {mustBeNumeric}
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

    % Convert Cartesian to spherical coordinates
    r = norm(p);
    theta = acos(p(3)/r);
    
    validatePosition(r, R);
    
    % Calculate wave properties
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
    % Compute the incident wave using a plane wave or point source
    if ~use_point_source
        P0 = exp(-1i*k0*r*mu);
    else
        % Monopoole source located in (0, 0, D)
        r0 = sqrt(D^2 + r^2 - 2*D*r*mu);
        P0 = exp(1i*k0*r0) / (1i * k0 * r0);
    end
end

function P = computeModalExpansion(order, k0, k1, R, rho0, c0, rho1, c1, kr, mu, use_point_source, D)
    % Compute the scattered field using modal expansion
    P = 0;
    for m = 0:order
        % Compute modal coefficient
        Am = computeModalCoefficient(m, k0, k1, R, rho0, c0, rho1, c1, use_point_source, D);
        
        % Compute spherical wave functions
        h1_m = computeSphericalHankel(m, kr);
        Pm = evaluateLegendrePolynomial(m, mu);
        
        % Add modal contribution
        P = P + Am * Pm * h1_m;
    end
end

function Am = computeModalCoefficient(m, k0, k1, r, rho0, c0, rho1, c1, use_point_source, D)
    % Compute the modal coefficient Am according to equation (14) in 
    % https://pmc.ncbi.nlm.nih.gov/articles/PMC3132099/
    
    if use_point_source
        Lm = computeSphericalHankel(m, k0*D);
    else
        Lm = (-1i)^m;
    end

    Z0 = rho0 * c0;
    Z1 = rho1 * c1;
    
    k0r = k0*r;
    k1r = k1*r;
    
    % Compute spherical Bessel functions and derivatives
    j1_m_k1r = computeSphericalBesseljDerivative(m, k1r);
    jm_k0r = computeSphericalBesselj(m, k0r);
    j1_m_k0r = computeSphericalBesseljDerivative(m, k0r);
    jm_k1r = computeSphericalBesselj(m, k1r);
    hm_k0r = computeSphericalHankel(m, k0r);
    h1_m_k0r = computeSphericalHankelDerivative(m, k0r);
    
    % Compute coefficient using Anderson's formula
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
    % Validate that the evaluation point is outside the sphere
    if r <= R
        error('AndersonSolution:InvalidPosition', ...
              'The evaluation point must be outside the sphere (r > R)');
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