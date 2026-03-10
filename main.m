%% MAIN
clc;
clear;
close all;
echo off; % Used to avoid unwanted warnings and other console stuff

addpath("functions\")

track = genTrack()
dist = getDistBounds(track.m,track)
len = getLineLength(track.m)

% disp(track.m)