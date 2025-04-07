function recon=Season_L2(data,flag,nd,lambda1,lambda2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reconstruction of Season_L2 framework
    % inputs:
    % data       raw data time-series
    % flag       quality flag for raw data(a flag value of 0 for good quality, 1 for uncertain quality and 3 for poor quality)
    % nd         one year length of time-series
    % lambda1    parameter of temporal neighborhood term
    % lambda2    parameter of interannual similarity term
    % outpus:
    % recon      reconstructed data time-series
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for more details,
    % see the following paper,
    % 储栋,管小彬,沈焕锋.考虑全时间序列信息的NDVI变分重建方法[J].地理与地理信息科学,2023,39(03):31-39.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~,~,DTD,LTL] = DAndLMatrixCompute(length(data),nd);
    w = zeros(size(flag));
    w(flag==0)=1;
    w(flag==1)=0.8;
    Weight = diag(w);
    equation_left=Weight+lambda1*DTD+lambda2*LTL;
    equation_right=Weight*data;
    [C,label]=chol(equation_left); %by Cholesky decomposition
    if ~label
        recon=C\(C'\equation_right);
    else
        recon=equation_left\equation_right;
    end
end
