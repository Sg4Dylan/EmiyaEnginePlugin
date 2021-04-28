classdef SignalUtils
    %SIGNALUTILS
    %   Audio processing
    
    properties
        EdgeFreq;
        ChanSum;
        AccCount;
    end
    
    properties (Constant)
        % EnvelopeDetect parameter
        TGain = 85;
        SlideWindow = 8;
        ChanSize = 1024;
    end
    
    methods (Access = private)
        function utils = EnvelopeDetect(utils,Input,InputSR,SrcFreq)
            % Must suit for fixed size
            if size(Input) ~= size(utils.ChanSum)
                return
            end
            % Prepare channel data
            if size(Input,2) == 2
                Chan = (Input(:,1) + Input(:,2))/2; % 转换到中置声道
            else
                Chan = Input(:,1);
            end
            % Update mean value
            NewChan = abs(fft(Chan));
            utils.ChanSum = (utils.AccCount.*utils.ChanSum+NewChan)/(utils.AccCount+1);
            % Calculate max value
            ChanMax = max(utils.ChanSum);
            % Find edge frequence
            ProbFreq = 0;
            for p = 1:utils.ChanSize
                SlideMean = 0;
                if (p+utils.SlideWindow) <= utils.ChanSize
                    SlideMean = mean(utils.ChanSum(p:p+utils.SlideWindow));
                else
                    SlideMean = mean(utils.ChanSum(p:utils.ChanSize))*(utils.SlideWindow/(utils.ChanSize-p));
                end
                if SlideMean < (ChanMax/pow2(utils.TGain/6))
                    ProbFreq = ProbFreq + p*InputSR/2/1024;
                    break
                end
            end
            % Check & fintune
            if utils.EdgeFreq < 1
                utils.EdgeFreq = round(SrcFreq);
            end
            if ProbFreq > 5000
                ProbFreq = ProbFreq-700;
            end
            if ProbFreq < 1
                return
            end
            ProbFreq = round(ProbFreq/2);
            utils.EdgeFreq = ProbFreq;
        end
        
        function Output = ModMonoSignal(~,Input,InputSR,SrcFreq,DstFreq,MGain)
            Output = Input;
        end
    end
    
    methods
        function utils = SignalUtils()
            utils.EdgeFreq = 0;
            utils.ChanSum = complex(zeros(utils.ChanSize,1));
            utils.AccCount = 0;
        end
        
        function Output = AkkoJitter(~,Input,JitLow,JitUp,DynPtc)
            % Generate random vector
            JitterVector = (JitUp-JitLow).*rand(size(Input))+JitLow;
            % Dynamic Protect
            if DynPtc == false
                Output = Input.*(1+JitterVector);
            else
                Adp = Input.*JitterVector;
                AdpPower = mean(abs(Adp),'all');
                SrcPower = mean(abs(Input),'all');
                SrcF = 1-(AdpPower/SrcPower);
                Output = Adp + SrcF.*Input;
            end
        end
        
        function Output = CopyBand(utils,Input,InputSR,SrcFreq,DstFreq,MGain)
            % Handle Mono/Stereo
            if size(Input,2) == 1
                Chan = Input(:,1);
            elseif size(Input,2) == 2
                Chan = (Input(:,1) + Input(:,2))/2; % Convert to Mid-Side
            else
                Output = Input;
                return
            end
            % Modulation
            ModOutput = utils.ModMonoSignal(Chan,InputSR,SrcFreq,DstFreq,MGain);
            % Return results
            if size(Input,2) == 1
                Output = ModOutput;
            else
                SideChan = Input(:,1) - Input(:,2);
                Output = [(ModOutput+SideChan)/2, (ModOutput-SideChan)/2];
            end
        end
        
        function utils = ResetEdgeFreq(utils)
            utils.EdgeFreq = 0;
            utils.ChanSum = complex(zeros(utils.ChanSize,1));
            utils.AccCount = 0;
        end
        
        function EdgeFreq = GetEdgeFreq(utils,Input,InputSR,SrcFreq)
            utils = EnvelopeDetect(utils,Input,InputSR,SrcFreq);
            EdgeFreq = utils.EdgeFreq;
        end
        
    end
end

