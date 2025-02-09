clc
clear 
close all

%load('data/TC12data.mat')
%writetable(data,'data/TC12data.csv')

load('data/TC13Data.mat')
writetable(data,'data/TC13data.csv')

load('data/allweather.mat')
writetable(data,'data/allweather.csv')