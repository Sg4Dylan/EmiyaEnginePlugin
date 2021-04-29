classdef SignalUtils < handle
    %SIGNALUTILS
    %   Audio processing
    
    properties (Access = private)
        EdgeFreq;
        ChanSum;
        AccCount;
    end
    
    properties (Constant)
        % EnvelopeDetect parameter
        TGain = 50;
        SlideWindow = 4;
        FFTSize = 512;
    end
    
    methods (Access = private)
        function utils = EnvelopeDetect(utils,Input,InputSR,DstFreq)
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
            utils.AccCount = utils.AccCount + 1;
            % Post process FFT result
            ChanFFT = utils.ChanSum/(2*utils.FFTSize);
            ChanFFT = 2*ChanFFT(1:utils.FFTSize);
            % Calculate max value
            ChanMax = max(ChanFFT);
            % Find edge frequence
            ProbFreq = 0;
            for p = 1:utils.FFTSize
                SlideMean = 0;
                if (p+utils.SlideWindow) <= utils.FFTSize
                    SlideMean = mean(ChanFFT(p:p+utils.SlideWindow));
                else
                    SlideMean = mean(ChanFFT(p:utils.FFTSize))*(utils.SlideWindow/(utils.FFTSize-p));
                end
                if SlideMean < (ChanMax/pow2(utils.TGain/6))
                    ProbFreq = (InputSR/2)/(p/utils.FFTSize);
                    break
                end
            end
            % Check & fintune
            if utils.EdgeFreq < 1
                utils.EdgeFreq = round(DstFreq);
            end
            if ProbFreq > 5000
                ProbFreq = ProbFreq-7000;
            end
            if ProbFreq < 1
                return
            end
            ProbFreq = round(ProbFreq/2);
            % Update EdgeFreq
            NewEdgeFreq = utils.EdgeFreq*0.25+ProbFreq*0.75;
            utils.EdgeFreq = NewEdgeFreq;
        end
        
        function Output = ModMonoSignal(~,Input,InputSR,SrcFreq,DstFreq)
            % Ref freq range
            FreqBegin = round(0.5+size(Input,1)*DstFreq/InputSR);
            FreqEnd = round(FreqBegin+0.0125*size(Input,1));
            % Src signal power
            SrcAbsSignal = abs(fft(Input));
            SrcPower = mean(SrcAbsSignal(FreqBegin:FreqEnd),'all');
            % Highpass filter
            [b,a] = butter(3,SrcFreq/(InputSR/2),'high');
            HighSignal = filter(b,a,Input);
            % Move Freq - Right
            SrcPos = round(SrcFreq*size(Input, 1)/InputSR);
            DstPos = round(DstFreq*size(Input, 1)/InputSR);
            EndPos = round(size(Input,1)/2);
            MoveSignal = fft(HighSignal);
            MoveSignal(SrcPos:SrcPos+EndPos-DstPos) = MoveSignal(DstPos:EndPos);
            MoveSignal(1:DstPos) = 0;
            % Move Freq -Left
            SrcPos = round(size(Input,1)-SrcPos);
            DstPos = round(size(Input,1)-DstPos);
            MoveSignal(EndPos:DstPos) = MoveSignal(SrcPos-DstPos+EndPos:SrcPos);
            MoveSignal(SrcPos:end) = 0;
            % Move power
            MoveAbsSignal = abs(MoveSignal);
            MovePower = mean(MoveAbsSignal(FreqBegin:FreqEnd),'all');
            % Calc gain
            Gain = SrcPower/MovePower;
            % Post highpass filter
            [b,a] = butter(5,DstFreq/(InputSR/2),'high');
            % Output results
            Output = filter(b,a,Gain.*real(ifft(MoveSignal)));
        end
    end
    
    methods
        function utils = SignalUtils()
            utils.EdgeFreq = 0;
            utils.ChanSum = complex(zeros(utils.FFTSize*2,1));
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
        
        function Output = CopyBand(utils,Input,InputSR,SrcFreq,DstFreq,DynPtc,OnlyModSignal)
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
            SrcF = 1.0;
            if DynPtc == true
                AdpPower = mean(abs(ModOutput),'all');
                SrcPower = mean(abs(MidChan),'all');
                if AdpPower/SrcPower < 1
                    SrcF = 1-(AdpPower/SrcPower);
                end
            end
            % Debug
            if OnlyModSignal
                MidChan = ModOutput;
            else
                MidChan = ModOutput + SrcF.*MidChan;
            end
            % Return results
            if size(Input,2) == 1
                Output = MidChan;
            else
                if OnlyModSignal
                    Output = [MidChan, MidChan];
                else
                    Output = [MidChan+SideChan, MidChan-SideChan];
                end
            end
        end
        
        function EdgeFreq = GetEdgeFreq(utils,Input,InputSR,DstFreq)
            utils.EnvelopeDetect(Input,InputSR,DstFreq);
            EdgeFreq = utils.EdgeFreq;
        end
    end
end

