function [yout, ystring] = poolData(yin,nVars,polyorder,usesine, laurentorder,  dyin, dyorder)
% Copyright 2017, All Rights Reserved
% Code by Niall Mangan for paper "Inferring biological networks by sparse
% identification of nonlinear dynamics"
% by N. M. Mangan S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ystring{ind} = '1';
ind = ind+1;



% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    yout(1,:);
    ystring{ind} = ['x_' num2str(i)];
    ind = ind+1;
    
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            yout(1,:);
            ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j) '}'];
            ind = ind+1;
            
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
            ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j) '}*x_{' num2str(k) '}'];
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                        '}*x_{' num2str(k) '}*x_{' num2str(l) '}'];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                        '}*x_{' num2str(k) '}*x_{' num2str(l)  '}*x_{' num2str(m) '}'];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(polyorder>=6)
    % poly order 6
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        for gg = m:nVars
                            yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m).*yin(:,gg);
                            ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                                '}*x_{' num2str(k) '}*x_{' num2str(l)  '}*x_{' num2str(m)...
                                '}*x_{' num2str(gg) '}'];
                            ind = ind+1;
                        end
                    end
                end
            end
        end
    end
end

if(polyorder>=7)
    % poly order 7
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        for gg = m:nVars
                            for pp = gg:nVars
                                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k)...
                                    .*yin(:,l).*yin(:,m).*yin(:,gg).*yin(:,pp);
                                ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                                    '}*x_{' num2str(k) '}*x_{' num2str(l)  '}*x_{' num2str(m)...
                                    '}*x_{' num2str(gg) '}*x_{' num2str(pp) '}'];
                                ind = ind+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

if(polyorder>=8)
    % poly order 8
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        for gg = m:nVars
                            for pp = gg:nVars
                                for qq = pp:nVars
                                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k)...
                                    .*yin(:,l).*yin(:,m).*yin(:,gg).*yin(:,pp)...
                                    .*yin(:,qq);
                                ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                                    '}*x_{' num2str(k) '}*x_{' num2str(l)  '}*x_{' num2str(m)...
                                    '}*x_{' num2str(gg) '}*x_{' num2str(pp)...
                                    '}*x_{' num2str(qq) '}'];
                                ind = ind+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if(polyorder>=9)
    % poly order 9
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        for gg = m:nVars
                            for pp = gg:nVars
                                for qq = pp:nVars
                                    for rr = qq:nVars
                                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k)...
                                    .*yin(:,l).*yin(:,m).*yin(:,gg).*yin(:,pp)...
                                    .*yin(:,qq).*yin(:,rr);
                                ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                                    '}*x_{' num2str(k) '}*x_{' num2str(l)  '}*x_{' num2str(m)...
                                    '}*x_{' num2str(gg) '}*x_{' num2str(pp)...
                                    '}*x_{' num2str(qq) '}*x_{' num2str(rr) '}'];
                                ind = ind+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if(polyorder>=10)
    % poly order 10
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        for gg = m:nVars
                            for pp = gg:nVars
                                for qq = pp:nVars
                                    for rr = qq:nVars
                                        for ss = rr:nVars
                                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k)...
                                    .*yin(:,l).*yin(:,m).*yin(:,gg).*yin(:,pp)...
                                    .*yin(:,qq).*yin(:,rr).*yin(:,ss);
                                ystring{ind} = ['x_{' num2str(i) '}*x_{' num2str(j)...
                                    '}*x_{' num2str(k) '}*x_{' num2str(l)  '}*x_{' num2str(m)...
                                    '}*x_{' num2str(gg) '}*x_{' num2str(pp)...
                                    '}*x_{' num2str(qq) '}*x_{' num2str(rr)...
                                    '}*x_{' num2str(ss) '}'];
                                ind = ind+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
    

thresh = 2e-7;

if(laurentorder>=1)
    disp('adding laurent order 1 terms')
    for i=1:nVars
        yout(yin(:,i)>=thresh,ind) = 1./yin(yin(:,i)>=thresh, i);
        yout(yin(:,i)<thresh, ind) = 0;
        ystring{ind} = ['1/x_' num2str(i)];
        ind = ind+1;
    end
end

if(laurentorder>=2)
    disp('adding laurent order 2 terms')
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            ijthresh = yin(:,i)>=thresh & yin(:,j)>=thresh;
            yout(ijthresh,ind) = 1./yin(ijthresh,i)./yin(ijthresh,j);
            yout(~ijthresh,ind) = 0;
            ystring{ind} = ['1/x_' num2str(i) '1/x_' num2str(j)];
            ind = ind+1;
            
        end
    end
end
if(laurentorder>=3)
    disp('adding laurent order 3 terms')
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                ijkthresh = yin(:,i)>=thresh & yin(:,j)>=thresh & yin(:,k)>=thresh;
                yout(ijkthresh,ind) = 1./yin(ijkthresh,i)./yin(ijkthresh,j)./yin(ijkthresh,k);
                yout(~ijkthresh, ind)= 0;
                ystring{ind} = ['1/x_{' num2str(i) '}1/x_{' num2str(j) '}1/x_{' num2str(k) '}'];
                ind = ind+1;
            end
        end
    end
end

if(usesine)
    for k=1:10;
        yout = [yout sin(k*yin) cos(k*yin)];
        warning('have not contsructed Thetastring for sin and cos')
    end
end

ynoder = yout;
size(yout);
% add term which are the first derivative (dy/dt)* all other terms tried.
if dyin~=0 & dyorder>0
    disp('yes adding derivatives!')
    yout;
    for jj = 1:size(dyin,2)
        for ii = 1:size(yout,2)
            dyout(:,ii) = dyin(:,jj).*yout(:,ii);
            ystring{ind} = [ystring{ii} 'dx_{' num2str(jj) '}/dt'];
            ind = ind+1;
        end
        yout = [yout dyout];
    end
end

% add square of derivatives (dy/dt)^2 * all other terms tried (avoid doing
% these twice as would happen if ynoder was replaced with yout calculated
% with first derivatives.

if dyin~=0 & dyorder >1
    for ii = 1:size(yout, 2)
        dyout(:,ii) = dyin.*dyin.*ynoder(:,ii);
    end
    yout = [yout dyout];
end

    