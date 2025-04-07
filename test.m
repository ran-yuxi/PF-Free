% Author: Yuxi ran
% Email: yx.ran@whu.edu.cn
% The codes is created based on the method described in the following paper
% A Parameter and Flag Adaptive Reconstruction Method for Satellite Vegetation Index Time Series.
%% a case of PF-Free method for reconstructing MOD13A2-NDVI time-series
%% I load data
load("NDVI.mat");
%% II flag self-generation 
nd = 23; % one year length of time-series
merge = 1; 
selfFlag = flagSelfGeneration(NDVI,nd,merge);
%% III parameters self-selection
[lambda,k] = parasSelfSelection(NDVI,selfFlag,nd);
%% IV reconstruction
reconNDVI = Season_L2(NDVI,selfFlag,nd,lambda,lambda*k);
%% V plot
figure("Color",'w');
DOY = 1:length(NDVI);
plot(DOY',NDVI,'b--');
hold on;
plot(DOY',reconNDVI,'m');
scatter(DOY(selfFlag==3),NDVI(selfFlag==3),'k','filled');
scatter(DOY(selfFlag==1),NDVI(selfFlag==1),'k');
scatter(DOY(selfFlag==0),NDVI(selfFlag==0),'r','filled');
hold off;
ylabel("NDVI");
legend("Original NDVI","reconNDVI","Good Points","Uncertain Points","Bad Points");


