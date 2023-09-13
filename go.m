tic
clc; clear; close all;
addpath('.\ExtractPlcr\')
addpath('.\PassiveAcousticPerception\')

% 动作和数字的数组
% actions = {'push', 'pull', 'left', 'right'};
% numbers = {'1', '2', '3', '4', '5'};

% 对所有的动作和数字进行排列组合
% for i = 1:length(actions)

%     for j = 1:length(numbers)
%         % 构造文件名
%         filename = [actions{i}, numbers{j}];
% 
%         % 执行你的函数
%         siz = getAcou(filename);
%         getPLCR(filename, siz);
%         fusion(filename);
%     end
% end
% toc

filename = 'N';
%siz = getAcou(filename);
getPLCR(filename, 7);
fusion(filename);
