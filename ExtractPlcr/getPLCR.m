function [] = getPLCR(FileName, siz)

time_diff1 = .5;
time_diff2 = 7; 
diff1 = (time_diff1 + time_diff2) ;
diff2 = (siz * 0.08) ;

DownSampleRate = 2;

% for trace_index = 1:length(TraceID)
% tempID = TraceID(trace_index);
% for dev_index = 1:length(device_Name)
dev_index = 1;
csi_data = [];
time_stamp = [];
csi = [];
csi_new = [];
tmp = 0;
%为方便一次提取多个，此处Filename为组合字符串，不需要的话自行改动

% if isempty(dir(['./Data/',FileName,'.mat'])) 

csi_trace = read_bf_file(['C:\Users\mengn\Desktop\summer\data\0721\',FileName,'.dat']);

% Initialize time_stamps
time_stamps = zeros(size(csi_trace));

% Iterate over csi_trace and convert timestamps to seconds.
for i = 1:10
    disp(csi_trace{i, 1});
    disp(class(csi_trace{i, 1}));  % Display the class of csi_trace{i}
    disp(isfield(csi_trace{i}, 'timestamp_low'));  % Check if 'timestamp_low' field exists
end
for i = 1:length(csi_trace)
    time_stamps(i) = csi_trace{i}.timestamp_low * 1e-6;
end

% Find the index of the first entry that is diff1 seconds after the start.
start_index = find(time_stamps >= time_stamps(1) + diff1, 1, 'first');

% Filter the csi_trace to remove first diff1 seconds.
filtered_csi_trace = csi_trace(start_index:end);

% Initialize filtered_time_stamps
filtered_time_stamps = zeros(size(filtered_csi_trace));

% Iterate over filtered_csi_trace and convert timestamps to seconds.
for i = 1:length(filtered_csi_trace)
    filtered_time_stamps(i) = filtered_csi_trace{i}.timestamp_low * 1e-6;
end

% Find the index of the last entry that is diff2 seconds after the start of filtered csi_trace.
end_index = find(filtered_time_stamps >= filtered_time_stamps(150) + diff2, 1, 'first');

% Filter the csi_trace to keep only first diff2 seconds.
final_csi_trace = filtered_csi_trace(1:end_index-1);
%final_time_stamps = filtered_time_stamps(1:end_index-1);
%final_csi_trace = filtered_csi_trace;
final_csi_trace = csi_trace;
csi_data = zeros(length(final_csi_trace), size(get_scaled_csi(final_csi_trace{1}), 2), size(get_scaled_csi(final_csi_trace{1}), 3));
time_stamp = zeros(length(final_csi_trace), 1);

for data_index = 1:length(final_csi_trace)
    csi_entry = final_csi_trace{data_index};
    csi_tmp = get_scaled_csi(csi_entry);
    csi_data(data_index, :, :) = squeeze(csi_tmp(1, :, :));
    time_stamp(data_index) = csi_entry.timestamp_low * 1e-6;
end

csi = csi_data(1:end,:,:);
csi_new(:,1:30) = csi_data(1:end,1,:);
csi_new(:,31:60) = csi_data(1:end,2,:);
csi_new(:,61:90) = csi_data(1:end,3,:);
csi_data = csi_new;
time_stamp = time_stamp(1:end)';
time_stamp = time_stamp - time_stamp(1,1);

%% 插值要对天线比插值 不能直接对原始的csi插值
all_time = ceil(100 * (time_stamp(1,end) - time_stamp(1,1)))/100;
NewTime = 0:0.001:all_time;

%save(['./MatData/',FileName],'csi_data','time_stamp','csi');
% else
%     load([FileName,'.mat']);
% end
csi_amplitude = mean(abs(csi));
csi_variance = sqrt(var(abs(csi)));

csi_ratio = squeeze(csi_amplitude ./ csi_variance);
csi_ratio = sum(csi_ratio,2); 
[~, midx] = max(csi_ratio);
% ant_ratio = mean(reshape(csi_ratio, F, A), 1);
%csi_value = squeeze(csi(:,:,:));
csi_value = squeeze(csi(:,:,:)./(csi(:,midx,:) + 0.01));
csi_value(:,midx,:) = [];



csi_value_down = csi_value(1:DownSampleRate:end,:,:);
time_stamp_down = time_stamp(:,1:DownSampleRate:end);

csi_value_new = interp1(time_stamp_down,csi_value_down,NewTime, 'spline');


csi_value_amp = csi_value_new;
[b_value,a_value] = DesignHPF(3,6,1000);
csi_value_amp = filter(b_value,a_value,csi_value_amp);
[b_value,a_value] = DesignLPF(70,80,1000);
csi_value_amp = filter(b_value,a_value,csi_value_amp);
csi_value = csi_value_amp;
for anti = 0:1
    csi_stft_new(:,anti * 30 + 1:anti * 30 + 30) = csi_value(:,anti + 1,:);
end
csi_stft_new(:,[1:5,26:35,56:60]) = [];

[com_value,f,t] = stft(csi_stft_new,1000,'Window',rectwin(201),'OverlapLength',101,'FFTLength',512);
com_value = sum(abs(com_value),3);
%com_value = abs(com_value(:,:,19));
%com_value = com_value(abs(f) < 50,:);
f_abs = abs(com_value);
%f_abs(4096,:) = 50;
[energyvalue,index] = max(f_abs,[],1);
plcr = - f(index) * 2.97e8/5.32e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plcr = hampel(plcr,5,1);
%plcr = hampel(plcr,5,0.5);
plcr = smooth(plcr,10);
plcr = sgolayfilt(plcr,3,13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengN(dev_index) = length(plcr);
DFS_Mat(1:length(plcr),dev_index) = plcr;
minLength = min(lengN);
DFS_Mat(minLength:end,:) = [];
lengN(dev_index) = length(plcr);
DFS_Mat(1:length(plcr),dev_index) = plcr;

DFS_Mat(1:length(plcr),dev_index+1) = t;
% end
minLength = min(lengN);
DFS_Mat(minLength:end,:) = [];
% figure;
% plot(DFS_Mat(:,1))
% ylim([-0.5,0.5]);

save(FileName,'DFS_Mat')
% end

savedData = DFS_Mat;
FileName = ['wifiRes', FileName, '.mat'];
save(FileName, 'savedData');




