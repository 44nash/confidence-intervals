

%%
%  Useful constants
%%
NUM_HIST_BINS= 40;

% size of each sample
SAMPLE_SIZE= 100;
%partOfSampleSize = 75;
%Zval =[.13, .25, .39, .52, .67, .84, 1.04, 1.28, 1.645];


%samPropMu = partOfSampleSize/SAMPLE_SIZE;
%disp(samPropMu)

%LowSideMU = samPropMu - Zval(1) * sqrt((samPropMu*(1-samPropMu))/SAMPLE_SIZE);
%HighSideMU = samPropMu + Zval(1) * sqrt((samPropMu*(1-samPropMu))/SAMPLE_SIZE); 

%
% Ground Truth 
%
% mu:  population mean
% 
% sigmaSq: population variance
%
MU= 57;
SIGMA= 15;
SIGMA_SQ= power(SIGMA,2);




MU_HAT_FIGURE=20;
SIGMA_SQ_FIGURE=21;


% the number of times sampled will be
% taken
TIME_INSTANTS= 400;

muHatValues = zeros(TIME_INSTANTS,1);
sigmaSqHatValues = zeros(TIME_INSTANTS,1);

centerMuHatValues = 0;
centerSigmaSqHatValues = 0;

%
%  Generate SAMPLE_SIZE many samples at time
%
%  t(1):  X(1,1), ..., X(1,SAMPLE_SIZE)
%
%  t(2):  X(2,1), ..., X(2,SAMPLE_SIZE)
%   .
%   .
%   .
%  t(m):  X(m,1), ..., X(m,SAMPLE_SIZE)
%
%  where m= TIME_INSTANTS
%
%  For each of m-many rounds, compute the
%  estimator muHat and sigmaSqHat, 
%  then display the sampling disribution

for t=1:1:TIME_INSTANTS,
  % take a smample ant time instant, t
  X = generateNormalSamples(MU,SIGMA_SQ, SAMPLE_SIZE);
  
  %compute MLE estimate for mu
  muHat = (1/SAMPLE_SIZE) * sum(X);
  
  %compute MLE estimate for sigmaSq
  sigmaSqHat = (1/SAMPLE_SIZE) * sum( power((X- muHat),2) );
%  disp("-----------------");
%  disp( sigmaSqHat);
  %record the MLE estimators
  muHatValues(t,1) = muHat;
  sigmaSqHatValues(t,1) = sigmaSqHat;
end

%disp("muhat Vals-------------- " +muHatValues);   ------------------------

muTText = 'sampling dist for $$\hat{\mu}$$';
muXText = '$$\hat{\mu}$$ values';
muYText = 'frequency';

[ muHatCenter, muHatStdErr]= visualizeSamplingDistribution(muHatValues,muTText,muXText,muYText,MU_HAT_FIGURE);
msg = sprintf('muHatCenter= %f  muHatStdErr=%f',muHatCenter,muHatStdErr);
disp(msg);

sigmaSqTText= 'sampling dist for $$\hat{\sigma}^{2}$$';
sigmaSqXText= '$$\hat{\sigma}^{2}$$ values';
sigmaSqYText= 'frequency';

[sigmaSqHatCenter, sigmaSqHatStdErr]= visualizeSamplingDistribution(sigmaSqHatValues,sigmaSqTText,sigmaSqXText,sigmaSqYText,SIGMA_SQ_FIGURE);
msg = sprintf('sigmaSqHatCenter=%f  sigmaSqHatStdErr=%f',sigmaSqHatCenter,sigmaSqHatStdErr);
disp(msg);

%-----------------------------Implentation----------------------------------------------------------------
SAMPLE_SIZE= 100;
partOfSampleSize = 42;
Zval =[.13, .25, .39, .52, .67, .84, 1.04, 1.28, 1.645];

x = randn(100,1);
%samPropMu = partOfSampleSize/SAMPLE_SIZE;
samPropMu = MU;
disp("-----------------"+samPropMu+"-------------------")
%Var = (1/SAMPLE_SIZE) * sum( power((SAMPLE_SIZE- samPropMu),2) );
Var = SIGMA_SQ;
disp(samPropMu)
for i=1:1:length(Zval),
percent = i*10;
LowSideMU = samPropMu - Zval(i) * sqrt((samPropMu*(1-samPropMu))/SAMPLE_SIZE);
HighSideMU = samPropMu + Zval(i) * sqrt((samPropMu*(1-samPropMu))/SAMPLE_SIZE); 

disp(" ");
disp("Confidence Intervals For MU "+percent+"%");
disp("LowSide "+LowSideMU );
disp("HighSide "+HighSideMU );




NUM_HIST_BINS=100;

centerHatValue = mean(x);

stdErrHatValue = sqrt(var(x));

%figure(i);
%clf;
%hist(x,NUM_HIST_BINS);
%y = linspace(0,1);
%plot(y,x,'-o')
%hold on;
%plot(samPropMu,0,'r^','MarkerSize',14);

%plot(LowSideMU,0,'g^','MarkerSize',14);
%plot(HighSideMU,0,'g^','MarkerSize',14);

%title("Confidence Intervals For	"+percent+"%");
%xlabel(xLabelText,'Interpreter','LaTex');
%ylabel(yLabelText);



end


disp("-------------------"+Var+"-------------------")

for i=1:1:length(Zval),
percent = i*10;
LowSideVar = Var - Zval(i) * sqrt((Var*(1-Var))/SAMPLE_SIZE);
HighSideVar = Var + Zval(i) * sqrt((Var*(1-Var))/SAMPLE_SIZE); 

disp(" ");
disp("Confidence Intervals For Var "+percent+"%");
disp("LowSide Var "+LowSideVar );
disp("HighSide Var "+HighSideVar );




NUM_HIST_BINS=100;

centerHatValue = mean(x);

stdErrHatValue = sqrt(var(x));

%figure(i+30);
%clf;
%hist(x,NUM_HIST_BINS);
%y = linspace(0,1);
%plot(y,x,'-o')
%hold on;
%plot(Var,0,'r^','MarkerSize',14);

%plot(LowSideVar,0,'g^','MarkerSize',14);
%plot(HighSideVar,0,'g^','MarkerSize',14);
 
%title("Confidence Intervals For Var	"+percent+"%");
%xlabel(xLabelText,'Interpreter','LaTex');
%ylabel(yLabelText);



end











