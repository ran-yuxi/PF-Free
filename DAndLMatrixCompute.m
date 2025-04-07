function [D,L,DTD,LTL] = DAndLMatrixCompute(n,nd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate the temporal neighborhood matrix D, interannual similarity matrix L, and D'*D, L'*L
    % inputs:
    % n       length of raw data time-series
    % nd      one year length of time-series
    % outputs:
    % D       temporal neighborhood matrix
    % L       interannual similarity matrix
    % DTD     DTD = D'*D
    % LTL     LTL = L'*L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for more details,
    % see the following papers:
    % 储栋,管小彬,沈焕锋.考虑全时间序列信息的NDVI变分重建方法[J].地理与地理信息科学,2023,39(03):31-39.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ny = n/nd;
    I2 = speye(ny*nd);
    O2 = zeros(ny*nd,1);
    D = [I2 O2 O2]+[O2 -2*I2 O2]+[O2 O2 I2];
    D = D(1:ny*nd,2:ny*nd+1);
    D(1,1) = -1;
    D(ny*nd,ny*nd) = -1;
    D = full(D);
    L = zeros(ny*nd);
    for i = 1:ny*nd-nd
        L(i,i) = 1;
        L(i,i+nd) = -1;
    end
    DTD = D'*D;
    LTL = L'*L;
end