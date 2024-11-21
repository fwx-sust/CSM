% Constants
EI = 1;  % Beam stiffness (arbitrary value)
c = 1;   % Constant load
x = linspace(0, 1, 100);
Q=0;M=0;

%shape functions
N1 = (1 - 3*x.^2 + 2*x.^3);
N2 = x .* (1 - x).^2;
N3 = (3*x.^2 - 2*x.^3);
N4 = x.^2 .* (1 - x);

%F
F = [0.5*c+Q; (1/12)*c-M; 0.5*c; (-1/12)*c];

%K
K = EI * [
    12,6,-12,6;
    6,4,-6,2;
   -12,-6,12,-6;
    6,2,-6,4];

% Apply boundary conditions (u(0) = 0, u'(0) = 0)
K_reduced = K(3:4, 3:4); % Remove rows/columns for fixed DOFs
F_reduced = F(3:4);      % Remove corresponding load vector entries

% Solve the reduced system
a_reduced = K_reduced \ F_reduced;

% Back-substitute constrained DOFs
a = zeros(4, 1);
a(3:4) = a_reduced;

% Construct the solution u_h(x)
u_h = N1 * a(1) + N2 * a(2) + N3 * a(3) + N4 * a(4);

% Plot the solution
figure;
plot(x, u_h)