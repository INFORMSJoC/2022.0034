function [ x_tilde ] = FarthestOnEllipsoid(origin, objective, A, b, c)

% Assumes ellipsoid defined by 1/2 * x'*A*x + b'*x + c <= 0

n = size(origin, 1);

%% Hessian Norm Maximization

% Second-order Taylor maximization
D = - objective.hess(origin);
e = objective.hess(origin)' * origin - objective.grad(origin);

L = chol(A)';

L_inv = inv(L);
C = L_inv * D * L_inv';

[Q, delta] = eig(C);

if any(any(imag(Q)))
    error('Imaginary eigenvalues')
else
    Q = real(Q);
end

if any(any(imag(delta)))
    error('Imaginary eigenvalues')
else
    delta = real(delta);
end

S = L_inv'* Q;
delta = diag(delta);

epsilon = S' * e;

alpha = ones(n, 1);

beta = S' * b;




%% Solve Hidden Convexity Problem


x = sdpvar(n, 1);
y = sdpvar(n, 1);

result = optimize([alpha'*y + beta'*x + c <= 0;
    0.5 * x.^2 - y <= 0], delta'*y + epsilon'*x, sdpsettings('solver', 'mosek', 'verbose', 0));

x_tilde = zeros(n, 1);
for i = 1 : n
    if abs( 0.5 * double(x(i))^2 - double(y(i)) ) < 1e-6
        x_tilde(i) = double(x(i));
    else
        x_tilde(i) = sqrt(2 * double(y(i)));
    end
end

x_tilde = S * x_tilde;

yalmip('clear')

end

