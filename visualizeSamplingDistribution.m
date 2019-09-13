function [centerHatValue, stdErrHatValue] = visualizeSamplingDistribution(theValues, titleText, xLabelText, yLabelText,figNum)

NUM_HIST_BINS=40;
%
% compute the center of the sampling dist for values
% 
centerHatValue = mean(theValues);

%
% compute the standard error of the sampling dist for muHat
%
stdErrHatValue = sqrt(var(theValues));

%
% visualize the sampling distributions for muHat 
% using histogram (mark the center with red triangle)
%
figure(figNum);
clf;
hist(theValues,NUM_HIST_BINS);
hold on;
plot(centerHatValue,0,'r^','MarkerSize',14);

title(titleText,'Interpreter','LaTex');
xlabel(xLabelText,'Interpreter','LaTex');
ylabel(yLabelText);


end