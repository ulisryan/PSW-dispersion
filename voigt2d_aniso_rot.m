function Z = voigt2d_aniso_rot(x, y, p)
%VOIGT2D_ANISO_ROT  Asymmetric, rotated 2D Voigt-like profile
%
%   Z = voigt2d_aniso_rot(x, y, p)
%
%   p = [z0, A, m, wxL, wyL, wxG, wyG, xc, yc, theta]
%       z0   : offset
%       A    : amplitude
%       m    : mixing (0 = pure Gaussian, 1 = pure Lorentzian)
%       wxL  : Lorentzian width in x'
%       wyL  : Lorentzian width in y'
%       wxG  : Gaussian  width in x'
%       wyG  : Gaussian  width in y'
%       xc,yc: center position
%       theta: rotation angle [rad]

z0   = p(1);
A    = p(2);
m    = p(3);
wxL  = p(4);
wyL  = p(5);
wxG  = p(6);
wyG  = p(7);
xc   = p(8);
yc   = p(9);
theta= p(10);

% shift to center
X = x - xc;
Y = y - yc;

% rotation
ct = cos(theta);
st = sin(theta);
xR =  ct.*X + st.*Y;
yR = -st.*X + ct.*Y;

% elliptic "radius" for L and G
rL2 = (xR./wxL).^2 + (yR./wyL).^2;
rG2 = (xR./wxG).^2 + (yR./wyG).^2;

% Lorentzian & Gaussian parts
L = 1 ./ (1 + 4*rL2);
G = exp(-4*log(2)*rG2);

% mixed Voigt-like profile
Z = z0 + A .* ( m.*L + (1-m).*G );
end
