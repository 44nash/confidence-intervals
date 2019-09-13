function [ X ] = generateNormalSamples(mu, sigmaSq, numSamples)
%
%  generateNormalSamples
%
%  produces samples from the Normal Distribution governed by
%  population parameters
%
%  mu- the mean (population parameter)
%
%  sigmaSq- the variance (population parameters)
%
%  numSamples- number of samples taken
%

X = normrnd(mu, sqrt(sigmaSq), numSamples,1);

end