function [Xi, indTheta, lambdavec, numterms, errorv] = ADMpareto(Theta, tol, pflag)

% Theta, a tolerance, and a plot flag
% calculates the nullspace of Theta with noise
% uses Donoho optimal shrinkage code to find the correct threshold for the
% singular values in the presence of noise.
% uses ADM algorithm to compute the linear combination of basis vectors
% spanning the nullspace of Theta that creates the sparsest resulting
% vector. This sparse vector gives the coefficients Xi so that Theta*Xi = 0
% and therefore select the terms in a sparse model.
% if the plot flag pflag =1 it plots the pareto front
% if pflag =2 then it will plot some diagnostic plots including the number
% of terms for each attempted initial condition. If you have some initial
% conditions that result in fairly sparse vectors, you can decrease the
% tolerance to improve the "resolution" of these from next best initial
% condtions
% returns 
% 1) Xi is a vector of the coefficients for the nonzero terms at each lambda
% 2) indTheta is a cell containing the indicies for the nonzero terms in
% Theta for ach value of lambda tried.
% 3)lambdavec, a lambda vector with the number of lambda values tried for that 
% variable's set of data 
% 4)numterms, a vector of the number of terms for each lambda value tried

Xi = zeros(size(Theta,2),1);
indTheta = cell(1,1);


%initial lambda value, which is the value used for soft thresholding in ADM
lambda = 1e-8;

jj = 1; % counter
num= 1; % initialize the number of nonzero terms found for the lambda

MaxIter = 1e4;

% use Donoho optimal shrinkage code to find null space in presence of
% noise.
[U, S, V]= svd(Theta', 'econ');
[m n] = size(Theta');
% m/n aspect ratio of matrix to be denoised
ydi = diag(Theta');
ydi(ydi< (optimal_SVHT_coef(m/n, 0)*median(ydi)))  = 0 ;
Theta2 = (U*diag(ydi)*V')';
nT = null(Theta2);


% vary lambda by factor of 2 until we hit the point
% where all coefficients are forced to zero.
% could also use bisection as commented out below
while num>0
     jj
    % use ADM algorithm with varying intial conditions to find coefficient
    % matrix
    [indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, pflag);
    
    indTheta{jj,1} = indTheta1; % save the indices of non zero coefficients
    Xi(:,jj) =Xi1; % get those coefficients
    numterms(jj,1) = numterms1;  % calculate how many terms are in the sparsest vector
  
    % calculate the error for the sparse vector found given this lambda
    errorv(jj,1) = sum((Theta*Xi(:,jj)));
    % store
    lambdavec(jj,1) = lambda;
    % index
    lambda = 2*lambda;
    num = numterms(jj,1);
    jj = jj+1;
    
    
if pflag>0
    figure(33)
    semilogy(numterms, abs(errorv),'o')
    hold on
    drawnow
    xlabel('Number of terms')
    ylabel('Error')
    title('Pareto Front')
    
    figure(34)
    loglog(lambdavec, numterms, 'o')
    hold on 
    drawnow
    xlabel('Lambda values')
    ylabel('Number of terms')
    
end
end


% 
% % use bisection method to find with higher
% % resolution where the drop to zero terms occured
%     lambda1= lambdavec(end-1);
%     lambda2 = lambdavec(end);
%     
%  kk = 1;
%  for kk = 1:3
%     
%     lambda_half = 0.5*(lambda2-lambda1);
%         % use ADM algorithm with varying intial conditions to find coefficient
%     % matrix
%     [indTheta1, Xi1, numterms1] = AMDinitvary(nT,lambda_half,MaxIter,tol, pflag);
%     
%     indTheta{jj,1} = indTheta1; % save the indices of non zero coefficients
%     Xi{jj,1} =Xi1; % get those coefficients
%     numterms(jj,1) = numterms1;  % calculate how many terms are in the sparsest vector
%     
%     % calculate the error for the sparse vector found given this lambda
%     errorv(jj,1) = sum(Theta(:, indTheta{jj,1})*Xi{jj,1});
%     % store
%     
%     % check which side the lambda_half fell on
%     lambdavec(jj,1) = lambda_half;
%     if numterms(jj,1) ==0
%         lambda2 = lambda_half;
%     else
%         lambda1 = lambda_half;
%     end
%     jj = jj+1;
%     %index
%  end
%  
%      
%  % now calculate Xi locally within the relavant order of
%  % magnitude of lambda
%  lambdafine = linspace(lambdavec(end)*0.1, lambdavec{ll,1}(end), 10);
% 
%  
 
 
        
    
