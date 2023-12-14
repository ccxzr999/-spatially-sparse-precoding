classdef channel_generation < handle
    %环境和优化类型的模拟
    properties
        numberTransmitAntennas
        numberRecieveAntennas
        numberDataStreams
        numberRFchains
        numberCluster
        numberRayPercluster
        angularSpread
    end
    properties (Access = public)
        channelMatrix
        arrayResponseTx
        arrayResponseRx
        alpha
    end
    
    methods
        function obj = channel_generation (parameters)
            obj.numberTransmitAntennas = parameters("numberTransmitAntennas");
            obj.numberRecieveAntennas = parameters("numberRecieveAntennas");
            obj.numberDataStreams = parameters("numberDataStreams");
            obj.numberRFchains = parameters("numberRFchains");
            obj.numberCluster = parameters("numberCluster");
            obj.numberRayPercluster = parameters("numberRayPercluster");
            obj.angularSpread = parameters("angularSpread");
            obj.launch ();
        end
        
        %%
        function obj = launch(obj)
            minAOD =-30;
            maxAOD =30;
            minEOD =-10;
            maxEOD =10;
            %簇服从均匀分布
            %transmit
            clusterAOD = rand(1,obj.numberCluster)*(maxAOD-minAOD)+minAOD;%[-30,30]
            clusterEOD = rand(1,obj.numberCluster)*(maxEOD-minEOD)+minEOD;%[-10,10]
            %reciever
            clusterAOA = (rand(1,obj.numberCluster)-0.5)*2*pi;%[-180,180]
            clusterEOA = (rand(1,obj.numberCluster)-0.5)*2*pi;%[-180,180]
            %径的分布拉普拉斯分布
            for k = 1:obj.numberCluster
                rayAOD((k-1)*obj.numberRayPercluster+1:k*obj.numberRayPercluster) = ...
                    laprnd(clusterAOD(k),obj.angularSpread,1,obj.numberRayPercluster);%发射端azimuth
                rayEOD((k-1)*obj.numberRayPercluster+1:k*obj.numberRayPercluster) = ...
                    laprnd(clusterEOD(k),obj.angularSpread,1,obj.numberRayPercluster);%发射端elevation
                
                rayAOA((k-1)*obj.numberRayPercluster+1:k*obj.numberRayPercluster) = ...
                    laprnd(clusterAOA(k),obj.angularSpread,1,obj.numberRayPercluster);%接收端azimuth
                rayEOA((k-1)*obj.numberRayPercluster+1:k*obj.numberRayPercluster) = ...
                    laprnd(clusterEOA(k),obj.angularSpread,1,obj.numberRayPercluster);%接收端elevation
            end
           %生成阵列响应
           obj.arrayResponseTx = zeros(obj.numberTransmitAntennas,obj.numberCluster*obj.numberRayPercluster);
           obj.arrayResponseRx = zeros(obj.numberRecieveAntennas,obj.numberCluster*obj.numberRayPercluster);
           for k = 1:obj.numberCluster*obj.numberRayPercluster
               %发天线阵列响应
               AT_row = zeros(1,sqrt(obj.numberTransmitAntennas));
               AT_col = zeros(1,sqrt(obj.numberTransmitAntennas));
               
               for w_t = 0:sqrt(obj.numberTransmitAntennas)-1
                   AT_row(w_t+1) = exp(1i*pi*(w_t*sin(rayAOD(k))*sin(rayEOD(k))));
               end
               for h_t = 0:sqrt(obj.numberTransmitAntennas)-1
                    AT_col(h_t+1) = exp(1i*pi*(h_t*cos(rayEOD(k))));
               end
               obj.arrayResponseTx(:,k) = (kron(AT_row,AT_col)./sqrt(obj.numberTransmitAntennas)).';
               %收天线阵列响应
               AR_row = zeros(1,sqrt(obj.numberRecieveAntennas));
               AR_col = zeros(1,sqrt(obj.numberRecieveAntennas));
               
               for w_t = 0:sqrt(obj.numberRecieveAntennas)-1
                   AR_row(w_t+1) = exp(1i*pi*(w_t*sin(rayAOA(k))*sin(rayEOA(k))));
               end
               for h_t = 0:sqrt(obj.numberRecieveAntennas)-1
                    AR_col(h_t+1) = exp(1i*pi*(h_t*cos(rayEOA(k))));
               end
               obj.arrayResponseRx(:,k) = (kron(AR_row,AR_col)./sqrt(obj.numberRecieveAntennas)).';
           end
           
            obj.generateComplexGainPath();
            obj.generateChannelMatrix();
        end
               
    %%
            function phi_c = laprnd(miu,sigma,m,n)
                beta = 1/(1-exp(sqrt(2)*pi/sigma));
                u = rand(m,n)-0.5;
                phi_c = miu - sigma*sign(u).*log(1-2*abs(u));
                phi_c = phi_c./beta;
            end
            
            function obj = generateComplexGainPath(obj)
                obj.alpha = sqrt(1/2)*(randn(1,obj.numberRayPercluster*obj.numberCluster)+...
                    1i*randn(1,obj.numberRayPercluster*obj.numberCluster));
            end
            
             function obj = generateChannelMatrix(obj)
                obj.channelMatrix = sqrt((obj.numberTransmitAntennas*obj.numberRecieveAntennas)/...
                    (obj.numberRayPercluster*obj.numberCluster))*obj.arrayResponseRx*diag(obj.alpha)*...
                    obj.arrayResponseTx';
             end
    end
end
        
        
                