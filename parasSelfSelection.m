function [lambda,k] = parasSelfSelection(data,flag,nd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adaptive determination of reconstruction parameters lambda and k using improved generalised cross-validation
    % inputs:
    % data       raw data time-series
    % flag       quality flag for raw data(a flag value of 0 for good quality, 1 for uncertain quality and 3 for poor quality)
    % nd         one year length of time-series
    % outputs:
    % lambda     parameter of temporal neighborhood term
    % k          Parameter for measuring the strength of temporal neighborhood term and interannual similarity term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for more details,
    % see the following papers:
    % Craven P, Wahba G. Smoothing noisy data with spline functions: estimating the correct degree of smoothing by the method of generalized cross-validation[J]. Numerische mathematik, 1978, 31(4): 377-403.
    % Brezinski C, Redivo-Zaglia M, Rodriguez G, et al. Multi-parameter regularization techniques for ill-conditioned linear systems[J]. Numerische Mathematik, 2003, 94: 203-228.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ny = length(data)/nd;
    w = zeros(size(flag));
    w(flag==0)=1; % 1 for good quality data point 
    w(flag==1)=0.8; % 0.8 for uncertain quality data point, 0 for bad quality data point
    [D,L,DTD,LTL] = DAndLMatrixCompute(length(data),nd);
    k = kCompute(data,w,D,L,DTD);
    H = HCompute(k,DTD,LTL);
    I = diag(sqrt(w));
    I = I(w~=0,:);
    [U,sigmai] = svdGCV(I,H);
    rH = nnz(w);
    lambda = GCV(data,w,U,sigmai,rH);
end
%% 
% adaptive parameter k function
function k = kCompute(NDVI,w,D,L,DTD)
    reconNDVI = double(wittakerFilter(NDVI,diag(w),DTD,1)); % Whittaker smoother
    Dx = sparse(D)*(reconNDVI);
    Lx = sparse(L)*(reconNDVI);
    k = (norm(Dx)^2)/(norm(Lx)^2);
    k(isnan(k) | isinf(k))=1;
end
%% 
% Whittaker smoother function
function z=wittakerFilter(y,Weight,DTD,lambda0)
    [C,flag] = chol(Weight + lambda0 * DTD);
    if ~flag
        z = (C\(C'\(Weight*y)));
    else
        z = (Weight + lambda0 * DTD)\(Weight*y);
    end
end
%% 
% Integration of the temporal neighborhood matrix with the interannual similarity matrix into one regularization matrix H
function [H] = HCompute(k,DTD,LTL)
    H0 = DTD+k*LTL;
    small_num = eps('single')*1e2;
    H0 = H0 + small_num*eye(size(H0));
    H = chol(H0);
end
%% 
% svd decomposition for GCV
function [U,sigmai] = svdGCV(I0,H0)
    [U,sigmai,~] = svd(I0/H0,"econ","vector");
    sigmai = sigmai(end:-1:1);
    U = U(:,end:-1:1);
end
%% 
% GCV function
function [lambda] = GCV(y,w,Ui,si,rH)
    y(w==0)=[];
    w(w==0)=[];
    ciH = (Ui.')*(sqrt(diag(w))*y);
    options = optimset('Display','off','Tolx',1e-4,'MaxFunEvals',200,'MaxIter',200);
    lb = 1e-2;
    ub = 2e3; 
    [lambda,fval,exitflag,output] = fminbnd(@GCVscore,lb,ub,options);
    function GCVs = GCVscore(s0)
        Ni = sqrt(sum(( (s0*ciH)./((si).^2+s0) ).^2 ));
        Di = sum((s0)./((si).^2+s0) );
        V = Ni/Di;
        GCVOri = V^2;
        GCVs = GCVOri;
    end
end