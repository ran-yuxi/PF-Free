function flag = flagSelfGeneration(data, nd, merge)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % flag self-generation
    % inputs:
    % data     raw data time-series
    % nd       one year length of of time-series
    % merge    Merge parameter, which considers whether to merge data when rearranging 
    % outputs: 
    % flag     quality flag for raw data(a flag value of 0 for good quality, 1 for uncertain quality and 3 for bad quality)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for more details,
    % see the following paper:
    % A Parameter and Flag Adaptive Reconstruction Method for Satellite Vegetation Index Time Series.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get the length of the time series, and the number of years
    len = length(data);
    ny = len/nd;
    num = (1:len).';

    % parameters preparation
    dataReStdPara = 0.1; % data standard deviation parameter
    dataReMinPara = 0.2; % data minimum value parameter（for feature points selection）
    dataReDisMinPara = 0.1; % minumum distance parameter
    pointsNumPara = 4; % available data point parameter
    dataLLPara = 0.1; % data minimum value parameter（for curve fitting）
    ratioParas = [0.6;0.4;0.2]; % ratio parameters
    fitNumberParas = [8;6;4;2]; % number of fitting parameters
    CTFParas = [0.025,0.05]; % range parameters
    CIParas = [0.25,0.05]; % confidence interval parameters

    % data rearrangement
    [dataRe, nd, ny] = rearrange(data, nd, ny, merge);
    dataRe2d = reshape(dataRe(1:(nd-1)*ny), ny, nd-1);
    dataReEnd = dataRe((nd-1)*ny+1:end);
    dataReMean = [mean(dataRe2d),mean(dataReEnd)];
    dataReStd = [std(dataRe2d),std(dataReEnd)];
    
    featurePointsRe2d = zeros(size(dataRe2d),'int8');
    featurePointsReEnd = zeros(size(dataReEnd),'int8');

    % feature points selection and process
    for od = 1:nd-1
        featurePoints = featurePointsSelection(dataRe2d(:,od),dataReStdPara,dataReMinPara,dataReDisMinPara,pointsNumPara);
        featurePoints = featurePointsProcess(dataRe2d(:,od),featurePoints,dataReStdPara);
        featurePointsRe2d(:,od) = featurePoints;
        clear featurePoints;
    end
    featurePointsReEnd = featurePointsSelection(dataReEnd,dataReStdPara,dataReMinPara,dataReDisMinPara,pointsNumPara);
    featurePointsReEnd = featurePointsProcess(dataReEnd,featurePointsReEnd,dataReStdPara);

    featurePointsRe = [reshape(featurePointsRe2d,size(featurePointsRe2d,1)*size(featurePointsRe2d,2),1); ...
          featurePointsReEnd];

    % check the result of feature points selection
    featurePointsPlot(dataRe,featurePointsRe,nd,ny,num);

    % curve fitting based feature points
    [fitParas, fitError, mu] = multiPolyFit(dataRe,featurePointsRe,nd,ny,num,ratioParas,fitNumberParas);
    [headCell, endCell] = bothEndConstantTermFit(dataRe,featurePointsRe,nd,ny,num,dataReMinPara,CTFParas);

    % create confidence interval
    [fitData, fitDelta] = CICreate(num,fitParas,fitError,mu,nd,ny,headCell,endCell,CIParas);

    % generate flag
    [flagRe, dataLLUL] = flagsGeneration(dataRe,fitData,fitDelta,num,dataLLPara);

    % check the result of curve fitting
    curveFittingPlot(dataRe,flagRe,fitData,dataLLUL,nd,ny,num);

    % flag inverse rearrangement
    [flag, nd, ny] = invRearrange(flagRe, nd, ny, merge);

    % check the result of generated flag
    flagPlot(data,flag,nd,ny,num);
end
%%
% data rearrangement function
% arranges the data into a single-year time series stacked over multiple years with the same period of data
function [dataRearrangement,ndNew,nyNew] = rearrange(data, nd, ny, merge)
    dataReshape = reshape(data,nd,ny);
    dataReshapeT = dataReshape';
    if (nargin<4 | merge==1) 
        merge = 1;
        ndNew = nd;
        nyNew = ny;
        dataRearrangement = reshape(dataReshapeT,nd*ny,1);
    else 
        ndNew = ceil(nd/merge);
        nyNew = ny*merge;
        dataRearrangement = zeros(nd*ny,1);
        
        for od = 1:ndNew-1
            dataReshapeTOd = dataReshapeT(:,merge*(od-1)+1:merge*od);
            dataRearrangement(nyNew*(od-1)+1:nyNew*od) = reshape(dataReshapeTOd',nyNew,1);
            clear dataReshapeTOd;
        end
        
        dataReshapeTOd = dataReshapeT(:,merge*(ndNew-1)+1:end); 
        dataRearrangement(nyNew*(ndNew-1)+1:end) = reshape(dataReshapeTOd',size(dataReshapeTOd,1)*size(dataReshapeTOd,2),1);
    end
end
%%
% feature points selection function
% filter feature points from each period of data
function featurePoints = featurePointsSelection(data,dataStdPara,dataMinPara,dataDisMinPara,pointsNumPara)
    dataMean = mean(data);
    dataStd = std(data);
    [dataSort,dataIndex] = sort(data,'descend'); 
    dataSortDiff = dataSort(1:end-1)-dataSort(2:end);
    dataSortDiffRatio = dataSortDiff./dataSort(1:end-1);
    dataSortDiffRatioStd = std(dataSortDiffRatio);
    
    featurePoints = zeros(size(data),'int8');
    if dataStd<dataStdPara 
        if dataMean<dataMinPara
            featurePointsIndex = [];
            featurePoints = zeros(size(data),'int8');
        else
            featurePointsIndex = linspace(1,length(data),length(data));
            featurePoints = ones(size(data),'int8');
        end
    else 
        dataSortDiffRatioFirst = dataSortDiffRatio(1);
        dataSortDiffRatioFirstIndex = 1;
        [dataSortDiffRatioHighest,dataSortDiffRatioHighestIndex] = max(dataSortDiffRatio);
        if (dataSortDiffRatioHighestIndex==1 | dataSortDiffRatioHighestIndex==2)
            featurePoints = zeros(size(data),'int8');
            featurePoints(dataIndex(1))=1; 
            featurePointsIndex = [1];
        else
            oriPoints = dataSortDiffRatio(dataSortDiffRatioFirstIndex:dataSortDiffRatioHighestIndex);
            oriPointsData = dataSort(dataSortDiffRatioFirstIndex:dataSortDiffRatioHighestIndex);
            oriPointsIndex = (dataSortDiffRatioFirstIndex:dataSortDiffRatioHighestIndex)';
            oriPointsLen = length(oriPointsIndex);
            dAbove = (oriPointsIndex-oriPointsLen).*(dataSortDiffRatioHighest-dataSortDiffRatioFirst)+(oriPoints-dataSortDiffRatioHighest)*(1-oriPointsLen);
            dBelow = sqrt((dataSortDiffRatioHighest-dataSortDiffRatioFirst).^2+(1-oriPointsLen).^2);
            d = abs(dAbove/dBelow);
            dGt = d(oriPointsData>dataDisMinPara);
            oriPointsGt = oriPoints(oriPointsData>dataDisMinPara);
            oriPointsDataGt = oriPointsData(oriPointsData>dataDisMinPara);
            oriPointsIndexGt = oriPointsIndex(oriPointsData>dataDisMinPara);
            if isempty(dGt) 
                featurePointsIndex = [];
                featurePoints = zeros(size(data),'int8');
            else
                if length(dGt)<=pointsNumPara 
                    [~,oriPointsDataGtMaxIndex] = max(oriPointsDataGt);
                    featurePoints = zeros(size(data),'int8');
                    dataSortDiffRatioEndIndex = oriPointsIndexGt(oriPointsDataGtMaxIndex);
                    featurePointsIndex = dataSortDiffRatioFirstIndex:dataSortDiffRatioEndIndex;
                    featurePoints(dataIndex(featurePointsIndex))=1; 
                else
                    [dGtSort,dGtIndex] = sort(dGt,'descend');
                    oriPointsDataGtSort = oriPointsDataGt(dGtIndex);
                    if dataSortDiffRatioStd <= dataStdPara
                        oriPointsDataGtSort = oriPointsDataGtSort(1:2);
                        dGtIndex = dGtIndex(1:2);
                    else
                        oriPointsDataGtSort = oriPointsDataGtSort(1:4);
                        dGtIndex = dGtIndex(1:4);
                    end
                    oriPointsDataGtMaxIndex = dGtIndex(find(oriPointsDataGtSort==max(oriPointsDataGtSort)));
                    dataSortDiffRatioEndIndex = oriPointsIndexGt(oriPointsDataGtMaxIndex);
                    featurePointsIndex = dataSortDiffRatioFirstIndex:dataSortDiffRatioEndIndex;
                    featurePoints(dataIndex(featurePointsIndex))=1; 
                end
            end
        end
    end
end
%%
% feature points processing function
% the standard deviation of the filtered feature points is greater than 0.1 
% and the feature points of the period that is not full of feature points are taken as the top 50 percent
function featurePointsAP = featurePointsProcess(data,featurePoints,dataStdPara)
    dataSelection = data(featurePoints==1);
    dataStd = std(dataSelection);
    if dataStd>dataStdPara & length(find(featurePoints==0))>=1;
        featurePointsSelection = find(featurePoints==1);
        featurePoints(featurePointsSelection(dataSelection<=prctile(dataSelection,50))) = 0;
    end
    featurePointsAP = featurePoints;
end
%%
% curve fitting function
% input feature points for polynomial curve fitting
function [fitParas,fitError,mu] = multiPolyFit(data,featurePoints,nd,ny,num,ratioParas,fitNumberParas)
    dataFit = data(featurePoints==1);
    numFit = num(featurePoints==1);
    if length(dataFit) > length(num)*ratioParas(1)
        [fitParas,fitError,mu] = polyfit(numFit,dataFit,fitNumberParas(1));
    else
        if length(dataFit) > length(num)*ratioParas(2)
            [fitParas,fitError,mu] = polyfit(numFit,dataFit,fitNumberParas(2));
        else
            if length(dataFit) > length(num)*ratioParas(3)
                [fitParas,fitError,mu] = polyfit(numFit,dataFit,fitNumberParas(3));
            else
                [fitParas,fitError,mu] = polyfit(numFit,dataFit,fitNumberParas(4));
            end
        end
    end
end
%%
% constant term fitting detection function at both ends 
% determine whether there are good data points at both ends
% if not, a constant term is needed to fit first
function [headCell,endCell] = bothEndConstantTermFit(data,featurePoints,nd,ny,num,dataMinPara,CTFParas)
    dataFit = data(featurePoints==1);
    numFit = num(featurePoints==1);
    if min(numFit)>ny
        headId = ceil(min(numFit)/ny);
        headValue = repmat(max(mean(dataFit(numFit>(headId-1)*ny & numFit<=headId*ny)),dataMinPara),ny*(headId-1),1);
        headIndex = linspace(1,ny*(headId-1),ny*(headId-1));
        headDelta_1 = repmat(CTFParas(1),ny*(headId-1),1);
        headDelta_2 = repmat(CTFParas(2),ny*(headId-1),1);
    else
        headValue = []; headId = []; headDelta_1 = []; headDelta_2 = [];
    end
    if max(numFit)<=(nd-1)*ny
        endId = ceil(max(numFit)/ny);
        endValue = repmat(max(mean(datFit(numFit>(endId-1)*ny & numFit<=endId*ny)),dataMinPara),ny*(nd-endId),1);
        endIndex = linspace(endId*ny+1,nd*ny,ny*(nd-endId));
        endDelta_1 = repmat(CTFParas(1),ny*(nd-endId),1);
        endDelta_2 = repmat(CTFParas(2),ny*(nd-endId),1);
    else
        endValue = []; endId = []; endDelta_1 = []; endDelta_2 = [];
    end
    headCell = {headValue,headId,headDelta_1,headDelta_2};
    endCell = {endValue,endId,endDelta_1,endDelta_2};
end
%%
% Confidence interval generator function 
% generate confidence intervals for the fitted curve
function [fitData,fitDelta] = CICreate(num,fitParas,fitError,mu,nd,ny,headCell,endCell,CIParas)
    [fitData_1,fitDelta_1] = polyconf(fitParas,num,fitError,'alpha',CIParas(1),'mu',mu);
    [fitData_2,fitDelta_2] = polyconf(fitParas,num,fitError,'alpha',CIParas(2),'mu',mu);
    if ~isempty(headCell{1})
        fitData_1(1:ny*(headCell{2}-1)) = headCell{1};
        fitData_2(1:ny*(headCell{2}-1)) = headCell{1};
        fitDelta_1(1:ny*(headCell{2}-1)) = headCell{3};
        fitDelta_2(1:ny*(headCell{2}-1)) = headCell{4};
    end
    if ~isempty(endCell{1})
        fitData_1(endCell{2}*ny+1:nd*ny) = endCell{1};
        fitData_2(endCell{2}*ny+1:nd*ny) = endCell{1};
        fitDelta_1(endCell{2}*ny+1:nd*ny) = endCell{3};
        fitDelta_2(endCell{2}*ny+1:nd*ny) = endCell{4};
    end
    fitData = {fitData_1,fitData_2};
    fitDelta = {fitDelta_1,fitDelta_2};
end
%%
% flag generation function 
% based on the spatial relationship between the data points and the confidence intervals of the fitted curves
function [flag,dataLLUL] = flagsGeneration(data,fitData,fitDelta,num,dataLLPara);
    dataLLParaList = dataLLPara.*ones(size(num));
    dataLL_1 = max(fitData{1}-fitDelta{1},dataLLParaList);
    dataUL_1 = min(max(fitData{1}+fitDelta{1},dataLLParaList),ones(size(num)));
    dataLL_2 = max(fitData{2}-fitDelta{2},dataLLParaList);
    dataUL_2 = min(max(fitData{2}+fitDelta{2},dataLLParaList),ones(size(num)));  
    
    dataBad = find(data<=dataLL_2 | data>dataUL_2);
    dataGood = find(data>dataLL_1 & data<=dataUL_1);
    dataUncer = find((data>dataUL_1 & data<=dataUL_2) | (data<=dataLL_1 & data>dataLL_2));

    flag = 3*ones(size(data));
    flag(dataBad)=3;
    flag(dataGood)=0;
    flag(dataUncer)=1;
    dataLLUL = {dataLL_1,dataUL_1,dataLL_2,dataUL_2};
end
%%
% inverse rearrangement function
function [data,ndNew,nyNew] = invRearrange(dataRe, nd, ny, merge)
    if (nargin<4 | merge==1)
        ndNew = nd;
        nyNew = ny;
        data = reshape(reshape(dataRe,ny,nd)',ny*nd,1);
    else 
        nyNew = ny/merge;
        ndNew = length(dataRe)/nyNew;
        endNum = (nd*merge-ndNew)*nyNew;
        data2d = zeros(ndNew,nyNew);
        for od = 1:nd-1
            data2d(merge*(od-1)+1:merge*od,:) = reshape(dataRe(ny*(od-1)+1:ny*od),merge,nyNew);
        end
        
        if endNum~=0 
            dataReEnd = dataRe(end-endNum+1:end);
            dataEnd = reshape(dataReEnd,length(dataReEnd)/nyNew,nyNew);
            data2d(merge*(nd-1)+1:end,:) = dataEnd;
        else
            data2d(merge*(nd-1)+1:end,:) = reshape(dataRe(ny*(nd-1)+1:ny*nd),merge,nyNew);
        end
        data = reshape(data2d,ndNew*nyNew,1);
    end
end
%%
% check the result of feature points selection
function featurePointsPlot(data,featurePoints,nd,ny,num)
    figure(1);
    plot(data);
    hold on;
    scatter(num(featurePoints==1),data(featurePoints==1),'r','filled');
    for i = 0:nd
        xline(ny*i,'k','LineWidth',0.8,'LineStyle','--');
    end
    hold off;
    legend("Original NDVI(rearrangeed)","feature Points");
end
%%
% check the result of curve fitting
function curveFittingPlot(data,flag,fitData,dataLLUL,nd,ny,num)
    figure(2);
    plot(data,'b');
    hold on;
    scatter(num(flag==0),data(flag==0),'r','filled');
    scatter(num(flag==1),data(flag==1),'k');
    scatter(num(flag==3),data(flag==3),'k','filled');
    plot(fitData{1},'m');
    plot(dataLLUL{1},'c--');
    plot(dataLLUL{2},'c--');
    plot(dataLLUL{3},'k--');
    plot(dataLLUL{4},'k--');
    for i = 0:nd
        xline(ny*i,'k','LineWidth',0.8,'LineStyle','--');
    end
    hold off;
    legend("Original NDVI(rearrangeed)","Good Points","Uncertain Points","Bad Points");
end
%%
% check the result of generated flag
function flagPlot(data,flag,nd,ny,num)
    figure(3);
    plot(data,'b');
    hold on;
    scatter(num(flag==0),data(flag==0),'r','filled');
    scatter(num(flag==1),data(flag==1),'k');
    scatter(num(flag==3),data(flag==3),'k','filled');
    for i = 0:nd
        xline(ny*i,'k','LineWidth',0.8,'LineStyle','--');
    end
    hold off;
    legend("Original NDVI","Good Points","Uncertain Points","Bad Points");
end