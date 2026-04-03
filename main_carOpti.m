%% MAIN
clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")

% Run Params
trackplot=true;

par = carParams();

track = genTrack(trackplot);
dist = getDistBounds(track.m,track);
len = getLineLength(track.m);

line = genLineTest(track.m);
time = calcLapTime(line);

