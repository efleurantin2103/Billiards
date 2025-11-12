function ax = evaluate_taylor(a, x)
% EVALUATETAYLOR - Evaluate tensor in Taylor polynomial basis
%
% SYNTAX:
%   ax = evaluatetaylor(a, x)
%
% INPUTS:
%   a : Tensor of Taylor coefficients
%       Dimensions: K_1 x ... x K_q x (N_1+1) x ... x (N_M+1)
%       where K_i are output dimensions and N_j is degree of variable j
%   x : M x Nx matrix, columns are evaluation points
%       M is number of Taylor variables
%
% OUTPUTS:
%   ax : K_1 x ... x K_q x Nx tensor of evaluated values
%
% DESCRIPTION:
%   Evaluates f(x) = sum a_alpha * x_1^alpha_1 * ... * x_M^alpha_M
%   where alpha_j ranges from 0 to N_j.

% Reshape x to 2D
sizex = size(x);
x = reshape(x, [sizex(1), prod(sizex(2:end))]);

% Handle intval types
if exist('intval', 'file')
    if isintval(a(1)) && ~isintval(x(1))
        x = intval(x);
    end
    if isintval(x(1)) && ~isintval(a(1))
        a = intval(a);
    end
end

MM = sizex(1);  % Number of Taylor variables

% Handle trivial case
if MM == 0
    ax = repmat(a, [1, sizex(2:end)]);
    return
end

% Determine tensor dimensions
da = length(size(a));
if da == 2 && size(a, 2) == 1
    da = 1;
end

% Handle scalar-valued case
if da == MM
    a = reshape(a, [1, size(a)]);
end

% Extract dimensions
sizea = size(a);
q = length(sizea) - MM;
Na = sizea(q+1:end);
N = Na - 1;

% Validate dimensions
if size(x, 1) ~= MM
    error('Dimension mismatch: size(x,1) must match number of variables');
end

Nx = size(x, 2);

% Create multi-index grid
index = cell(1, MM);
for m = 1:MM
    index{m} = 0:N(m);
end

[alpha{1:MM}] = ndgrid(index{:});
for m = 1:MM
    alpha{m} = alpha{m}(:);
end

% Evaluate Taylor monomials
xalpha = ones(prod(Na), 1);
for m = 1:MM
    am = length(alpha{m});
    xalpha = xalpha .* repmat(x(m, :), am, 1) .^ repmat(alpha{m}, 1, Nx);
end

% Reshape and multiply with coefficients
xalpha = reshape(xalpha, [ones([1, q]), Na, Nx]);
xalpha = repmat(xalpha, [sizea(1:q), ones(1, MM+1)]);

ax = a .* xalpha;
ax = reshape(ax, [sizea(1:q), prod(Na), Nx]);
ax = sum(ax, q+1);
ax = reshape(ax, [sizea(1:q), sizex(2:end)]);

% Determine tensor dimensions
da = length(size(a));
if da == 2 && size(a, 2) == 1
    da = 1;  % Column vector is 1D
end

% Handle scalar-valued case (q=1, K_1=1)
if da == MM
    a = reshape(a, [1, size(a)]);
end

% Trivial case: constant function
if MM == 0
    ax = repmat(a, [1, sizex(2:end)]);
    return
end

% Extract dimensions
sizea = size(a);
q = length(sizea) - MM;          % Number of output dimensions
K = sizea(1:q);                  % Output tensor dimensions
Na = sizea(q+1:end);             % Number of coefficients per variable
N = Na - 1;                      % Maximum degree for each variable

% Validate input dimensions
if size(x, 1) ~= MM
    error('Input dimension mismatch: size(x,1) must equal M');
end

Nx = size(x, 2);  % Number of evaluation points

% Create multi-index grid for Taylor basis
index = cell(1, MM);
for m = 1:MM
    index{m} = 0:N(m);  % Powers from 0 to maximum degree
end

[alpha{1:MM}] = ndgrid(index{:});
for m = 1:MM
    alpha{m} = alpha{m}(:);  % Flatten each index array
end

% Evaluate Taylor monomials: x_1^{alpha_1} * ... * x_M^{alpha_M}
xalpha = ones(prod(Na), 1);
for m = 1:MM
    am = length(alpha{m});
    % Compute x_m^{alpha_m} for all multi-indices and evaluation points
    xalpha = xalpha .* repmat(x(m, :), am, 1) .^ repmat(alpha{m}, 1, Nx);
end

% Reshape basis evaluations to match coefficient tensor
xalpha = reshape(xalpha, [ones([1, q]), Na, Nx]);
xalpha = repmat(xalpha, [sizea(1:q), ones(1, MM+1)]);

% Multiply coefficients by basis functions and sum
ax = a .* xalpha;
ax = reshape(ax, [sizea(1:q), prod(Na), Nx]);
ax = sum(ax, q+1);
ax = reshape(ax, [sizea(1:q), sizex(2:end)]);

end