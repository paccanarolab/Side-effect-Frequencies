function [ W, H, J ] = DecompositionAlgorithm( R, k, alpha )
% DECOMPOSITION ALGORITHM implements the method described in our
%paper: 
%   Galeano et al. Predicting the Frequency of Drug Side effects.
%   Input arguments:
%       * R: matrix of drugs x side effects with rating values.
%       * k: number of latent features.
%       * alpha: confidence on the zeros.
%       * maxiter: maximun number of iterations.
% Stopping criteria: until termination tolerance tolx or maxiter is
% satisfied.
% Copyright (C) 2017 by Diego Galeano.


    % Termination tolerance on change in the elements of W and H.
    tolx = 1e-3;   % default is 1e-3.
    maxiter = 2000; % default.
    
    % initial variance
    variance = 0.01;

    % Dimensions of the drug-side effect matrix
    [ndrugs, nses] = size(R);

    % Initialization
    W0 = rand(ndrugs, k)*sqrt(variance);
    H0 = rand(k, nses)*sqrt(variance);

    % normalization
    H0 = H0./repmat(sqrt(sum(H0.^2,2)),1,nses);

    % epsilon based on machine precision.
    sqrteps = sqrt(eps);

    % cost function values
    J = zeros(maxiter, 1);

    % filter for clinical trials values
    CT = R > 0;
    % filter for unobserved associations
    UN = R == 0;

    for iter = 1:maxiter

        numer = (CT.*R)*H0';

        W = max(0,W0 .* (numer./((CT .* (W0*H0) + alpha*UN .* (W0*H0))*H0' + eps(numer))));

        numer = W'*(CT.*R);
        H = max(0,H0.* (numer./ (W'*(CT.*(W*H0) + alpha*UN.*(W*H0))  + eps(numer))));


        % Compute cost function
        J(iter) =  0.5*norm(CT.*(R - W*H), 'fro')^2 +...
            0.5*alpha*(norm(UN.*(R - W*H), 'fro')^2);


        % Get norm of difference and max change in factors
        dw = max(max(abs(W-W0) / (sqrteps+max(max(abs(W0))))));
        dh = max(max(abs(H-H0) / (sqrteps+max(max(abs(H0))))));
        delta = max(dw, dh);

        % Check for convergence
        if iter > 1
            if delta <= tolx

                fprintf('\n iter %d delta %e\n', iter, delta);
                J(iter+1:end) = [];
                break;


            end
        end

        % Remember previous iteration results
        W0 = W;
        H0 = H./repmat(sqrt(sum(H.^2,2)),1,nses); % normalize
        
        fprintf('\n iter %d cost function %e\n', iter, J(iter));
    end

             
            

end
