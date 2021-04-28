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
                Chan = (Input(:,1) + Input(:,2))/2; % Convert to Mid-Side
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
        
        function Output = ModMonoSignal(~,Input,InputSR,SrcFreq,DstFreq)
            % Ref freq range
            FreqBegin = round(0.5+2*size(Input,1)*SrcFreq/InputSR);
            FreqEnd = round(FreqBegin+0.025*size(Input,1));
            % Src signal power
            SrcAbsSignal = abs(fft(Input));
            SrcPower = mean(SrcAbsSignal(FreqBegin:FreqEnd),'all');
            % Highpass filter
            [z,p] = butter(3,2*SrcFreq/InputSR,'high');
            HighSignal = filtfilt(z,p,Input);
            % Move freq
            MoveRange = round(size(Input,1)*(DstFreq-SrcFreq)/(InputSR/2));
            MoveSignal = circshift(fft(HighSignal),MoveRange);
            MoveSignal(1:round(size(Input,1)*DstFreq/(InputSR/2)-1)) = 0;
            % Move power
            MoveAbsSignal = abs(MoveSignal);
            MovePower = mean(MoveAbsSignal(FreqBegin:FreqEnd),'all');
            % Calc gain
            Gain = SrcPower/MovePower;
            % Output results
            Output = Gain.*real(ifft(MoveSignal));
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
        
        function Output = CopyBand(utils,Input,InputSR,SrcFreq,DstFreq,DynPtc)
            % Handle error
            if SrcFreq >= DstFreq || SrcFreq < 1 || DstFreq < 1 || SrcFreq > (InputSR/2) || DstFreq > (InputSR/2)
                Output = Input;
                return
            end
            % Handle Mono/Stereo
            if size(Input,2) == 1
                MidChan = Input(:,1);
            elseif size(Input,2) == 2
                MidChan = (Input(:,1) + Input(:,2))/2; % Convert to Mid-Side
                SideChan = (Input(:,1) - Input(:,2))/2;
            else
                Output = Input;
                return
            end
            % Modulation
            ModOutput = utils.ModMonoSignal(MidChan,InputSR,SrcFreq,DstFreq);
            % Dynamic Protect
            if DynPtc == true
                AdpPower = mean(abs(ModOutput),'all');
                SrcPower = mean(abs(MidChan),'all');
                SrcF = 1-(AdpPower/SrcPower);
                MidChan = ModOutput + SrcF.*MidChan;
            end
            % Return results
            if size(Input,2) == 1
                Output = MidChan;
            else
                Output = [MidChan+SideChan, MidChan-SideChan];
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

