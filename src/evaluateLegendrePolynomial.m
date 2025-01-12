function L = evaluateLegendrePolynomial(degree, x)
    % EVALUATELEGENDREPOLYNOMIAL Evaluates the Legendre polynomial of specified degree
    %   L = evaluateLegendrePolynomial(degree, x) evaluates the Legendre polynomial
    %   of degree n at point(s) x using the coefficient representation.
    %
    % Arguments:
    %   degree: Non-negative integer specifying the degree of the Legendre polynomial
    %   x: Scalar or vector of points at which to evaluate the polynomial
    %
    % Returns:
    %   L: Value(s) of the Legendre polynomial at specified point(s)
    %
    % Example:
    %   x = linspace(-1, 1, 100);
    %   L = evaluateLegendrePolynomial(2, x);
    %   plot(x, L);
    %   title('Legendre Polynomial P_2(x)');
    %
    % ABOUT:
    %     author       - Antonio Stanziola
    %     date         - 11th January 2025
    %     last update  - 11th January 2025
    
    validateattributes(degree, {'numeric'}, ...
        {'integer', 'nonnegative', 'scalar'}, ...
        'evaluateLegendrePolynomial', 'degree');
    
    coefficients = generateLegendrePolynomialCoefficients(degree);
    L = polyval(coefficients, x);
end

function coefficients = generateLegendrePolynomialCoefficients(n)
    % GENERATELEGENDREPOLYNOMIALCOEFFICIENTS Generate coefficients for Legendre polynomials
    %   coefficients = generateLegendrePolynomialCoefficients(n) generates the
    %   coefficients of the nth Legendre polynomial in descending powers of x.
    %
    % Arguments:
    %   n: Non-negative integer specifying the degree of the Legendre polynomial
    %
    % Returns:
    %   coefficients: Column vector of coefficients where the mth element is the
    %                 coefficient of x^(n+1-m)
    %
    % Algorithm:
    %   Uses the recurrence relation for Legendre polynomials:
    %   (k)P_k(x) = (2k-1)xP_{k-1}(x) - (k-1)P_{k-2}(x)

    validateattributes(n, {'numeric'}, ...
        {'integer', 'nonnegative', 'scalar'}, ...
        'generateLegendrePolynomialCoefficients', 'n');
    
    % Handle base cases
    if n == 0
        coefficients = 1;
        return
    elseif n == 1
        coefficients = [1; 0];
        return
    end
    
    % Initialize vectors for recurrence relation
    coefficients_k_minus_2 = zeros(n + 1, 1);  % P_{k-2}
    coefficients_k_minus_2(n + 1) = 1;
    
    coefficients_k_minus_1 = zeros(n + 1, 1);  % P_{k-1}
    coefficients_k_minus_1(n) = 1;
    
    % Apply recurrence relation
    for k = 2:n
        coefficients_k = zeros(n + 1, 1);
        
        % Calculate coefficients for current degree
        for exponent = n-k+1:2:n
            coefficients_k(exponent) = (2*k-1) * coefficients_k_minus_1(exponent+1) + ...
                                     (1-k) * coefficients_k_minus_2(exponent);
        end
        
        % Handle the constant term
        coefficients_k(n+1) = coefficients_k(n+1) + ...
                             (1-k) * coefficients_k_minus_2(n+1);
        
        % Normalize by k
        coefficients_k = coefficients_k / k;
        
        % Update coefficient vectors for next iteration
        if k < n
            coefficients_k_minus_2 = coefficients_k_minus_1;
            coefficients_k_minus_1 = coefficients_k;
        end
    end
    
    coefficients = coefficients_k;
end