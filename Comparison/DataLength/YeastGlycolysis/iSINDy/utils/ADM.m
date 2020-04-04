% Code from Finding a sparse vector in a subspace:
% Linear sparsity using alternating directions
% by Qing Qu, Ju Sun, and John Wright
% code found here  https://sites.google.com/site/homeqingqu/miscellaneous
% solve the following problem
% min_{q,x} 1/2*||Y*q - x||_2^2 + lambda * ||x||_1, s.t. ||q||_2 = 1
% by alternating minimization method (ADM).
% Y: input data
% q_init: initialization for q, lambda
% lambda: penalty parameter
% MaxIter: max iteration
% tol: tolerance for convergence
% q: output result
%% 
% Modified By: K.Kahirman
% Last Updated: 2019/07/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = ADM(Y,q_init,lambda,MaxIter,tol)
q = q_init;

for k = 1:MaxIter
    q_old = q;
    Dummy=Y*q;
    x=sign(Dummy).*max(abs(Dummy)-lambda,0);% update y by soft thresholding
    q = Y'*x/norm(Y'*x,2); % update q by projection to the sphere
    res_q = norm(q_old-q,2);
    
    if (res_q<=tol)
        return;
    end
    
    k;
end

