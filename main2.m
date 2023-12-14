clear all
%snrValuesDB = -40:5:0; %in Figures
snrValuesDB = 0; %in Figures
snrValues = 10.^(snrValuesDB./10);
tryNumber = 100; 
numberDataStreams = 1:1:4;
%% Parameters
parameters_main2 = containers.Map('KeyType','char','ValueType','any');%创建了一个新的containers.Map实例，并指定了键类型（KeyType）为字符（char），值类型（ValueType）为任意类型（any）。
parameters_main2("numberTransmitAntennas") = 64; % Number of transmit antennas parameters_main['key1'] = 'value1';  % 使用[]操作符插入键值对
parameters_main2("numberRecieveAntennas") = 16; % Number of receive antennas
parameters_main2("numberDataStreams") = 1; % Number of data streams
parameters_main2("numberRFchains") = 4; % Number of RF chains for precoding and combining
parameters_main2("numberCluster") = 8; % Number of clusters
parameters_main2("numberRayPercluster") = 10; % Number of rays per cluster
parameters_main2("angularSpread") = 7.5; % Angular spread of 7.5 degree
spectralEffOptimal = zeros(tryNumber,length(numberDataStreams));
spectralEffHybrid = zeros(tryNumber,length(numberDataStreams));
for s = 1 :length(numberDataStreams)
    SNR = snrValues;
    parameters_main2("numberDataStreams") = numberDataStreams(s); % Number of data streams
    for i = 1:tryNumber
        channel = channel_generation(parameters_main2);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,s) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,s) = tempObj.spectralEfficiency;
    end
end

spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
% figure(); 
hold on
l1 = plot(numberDataStreams,spectralEffOptimalSNR,'-s','color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);
l2 = plot(numberDataStreams,spectralEffHybridSNR,'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
legend([l1,l2],'Optimal unconstrained precoding N_s','Hybrid precoding and combining N_s','Location','northwest','FontSize',15);
xlabel("SNR (dB)",'FontSize', 20)
ylabel("Spectral Efficiency(bits/s/Hz)",'FontSize', 20)