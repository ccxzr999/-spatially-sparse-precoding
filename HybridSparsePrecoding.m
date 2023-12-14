classdef HybridSparsePrecoding < handle
    
    properties
        channelSimulation
        SNR
        spectralEfficiency
    end
    
    methods
        function obj = HybridSparsePrecoding(SNR,Channel)
            obj.channelSimulation = Channel;
            obj.SNR = SNR;
            obj.launch();
        end
        
        function obj = launch(obj)
            %%%% Algorithm 1 - Hybrid precoding
            H = obj.channelSimulation.channelMatrix;
            [~,~,V] = svd(H);
            fOpt = V(:,1:obj.channelSimulation.numberDataStreams);
            fRF = [];
            fRes = fOpt;
            for r = 1:obj.channelSimulation.numberRFchains
                psi = obj.channelSimulation.arrayResponseTx'*fRes;
                [~,k] = max(diag(psi*psi'));
                fRF = [fRF obj.channelSimulation.arrayResponseTx(:,k)];
                F_BB = (fRF'*fRF)\(fRF'*fOpt);
                fRes = (fOpt-fRF*F_BB)/norm(fOpt-fRF*F_BB,'fro');
            end
            F_BB = sqrt(obj.channelSimulation.numberDataStreams)*(F_BB/norm(fRF*F_BB,'fro'));
            % Algorithm 2 - Hybrid combining
            CovRx = (obj.SNR/obj.channelSimulation.numberDataStreams)*...
                H*fRF*F_BB*F_BB'*fRF'*H'+eye(obj.channelSimulation.numberRecieveAntennas);
            W_MMSE = ((1/sqrt(obj.SNR))*(F_BB'*fRF'*H'*H*fRF*F_BB+...
                (obj.channelSimulation.numberDataStreams/obj.SNR)*eye(obj.channelSimulation.numberDataStreams))...
                \ (F_BB'*fRF'*H'))';
             W_RF = [];
            W_Res = W_MMSE;
            for r = 1:obj.channelSimulation.numberRFchains
                psi = obj.channelSimulation.arrayResponseRx'*CovRx*W_Res;
                [~,k] = max(diag(psi*psi'));
                W_RF = [W_RF obj.channelSimulation.arrayResponseRx(:,k)];
                W_BB = (W_RF'*CovRx*W_RF)\(W_RF'*CovRx*W_MMSE);
                W_Res = (W_MMSE-W_RF*W_BB)/norm(W_MMSE-W_RF*W_BB,'fro');
            end
            Rn_Hybrid = W_BB'*W_RF'*W_RF*W_BB;
            obj.spectralEfficiency = abs(log2(det(eye(obj.channelSimulation.numberDataStreams)...
                + (obj.SNR/obj.channelSimulation.numberDataStreams)...
                * (Rn_Hybrid\(W_BB'*W_RF'*H*fRF*F_BB*F_BB'*fRF'*H'*W_RF*W_BB)))));
        end
    end
end
            
                
                