function [lambda] = scree_plot(Y)
% Scree plot

N = size(Y,2);
Y_bar = mean(Y,2);

% Empirical covariance
Yc = Y - Y_bar(:,ones(1,N));
S = Yc*(Yc')/N;

% Eigenvalues
lambda = svd(S);

% Scree plot
figure('Name','Scree plot','NumberTitle','Off');
title('Scree plot')
plot(log10(lambda));
ylabel('$\log_{10} (\lambda)$','Interpreter','LaTex','FontSize',12)
xlabel('$n$','Interpreter','LaTex','FontSize',12)

end
