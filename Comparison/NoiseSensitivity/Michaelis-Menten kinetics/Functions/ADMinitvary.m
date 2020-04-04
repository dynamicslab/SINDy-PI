function [indTheta, Xi, numterms] = ADMinitvary(nT,lambda,MaxIter,tol, plotflag)

% normalize the collumns of the null space of Theta to use in the
% initial conditoin for ADM search
for ii = 1:size(nT,2)
    nTn(:,ii) = nT(:,ii)/mean(nT(:,ii));
end

% run ADM algorithm on each row of nTn
for jj = 1:size(nTn,1)
    jj;
    qinit = nTn(jj,:)'; % intial conditions
    q(:,jj) = ADM(nT, qinit, lambda, MaxIter, tol); % algrorithm for
    % finding coefficients resutling in  sparsest vector in null space
    out(:,jj) = nT*q(:,jj); % compose sparsets vectors
    nzeros(jj)= length(find(abs(out(:,jj))<lambda)); % chech how many zeros each
    % of the found sparse vectors have
end

indsparse = find(nzeros==max(nzeros)); % find the vector with the largest number of zeros
indTheta = find(abs(out(:,indsparse(1)))>=lambda); % save the indices of non zero coefficients
% Xi = out(indTheta, indsparse(1)); % get those coefficients
Xi = out(:, indsparse(1)); % get sparsest vector
smallinds =(abs(out(:,indsparse(1)))<lambda);
Xi(smallinds)=0; % set thresholded coefficients to zero


% check that the solution found by ADM is unique.
if length(indsparse)>1
    Xidiff= (out(indTheta, indsparse(1))-out(indTheta,indsparse(2)));
    if Xidiff>tol
        warning('ADM has discovered two different sparsest vectors')
    end
end
% calculate how many terms are in the sparsest vector
numterms = length(indTheta);

if plotflag == 2
    figure(121)
    semilogy(size(nT,1)-nzeros, 'o');
    drawnow 
end

