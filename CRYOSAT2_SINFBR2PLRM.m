% CRYOSAT2_SINFBR2PLRM.M : is a MATLAB program that converts SARIN mode FBR
% data to PLRM (pseudo-LRM) waveforms. 
%
% Copyright (C) 2013 Natalia Galin.
% Questions/Comments/Suggestion/Corretions? 
% Please do get in touch: natalia.galin@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%-----------------------------------------------------------------------%%
%% INPUTS REQUIRED:
% filename - name of the CryoSat-2 SARIN mode FBR file you'd like to 
% process. 
% 
%% DEPENDENCIES:
%  Cryo_L1b_readv2.m - Version 2.0 modified to support SARIN FBR data. 
%  Function originally written by S. Dinardo (and
%  available for download here: https://earth.esa.int/web/guest/
%  software-tools/-/article/cryosat-matlab-routines). 
% 
%% OUTPUTS: 
%  PLRM(RECORD,BURST,SAMPLES) - the pulse limited waveforms from the SARIN
%  FBR mode data. (Note, these waveforms will be 512x2 samples long, not
%  the 128 samples as the native CryoSat-2 LRM waveforms.)
%
%  HDR - header information extracted from 'filename'.
%
%  CS  - data information extracted from 'filename'. Contains useful info
%  like latitude/longitude values to locate the PLRM waveforms.
%
%%-----------------------------------------------------------------------%% 

function [PLRM,HDR,CS] = CRYOSAT2_SINFBR2PLRM(filename)

[HDR, CS]=Cryo_L1b_readv2(filename);

CHANNEL1 = CS.FBR.Rx1;
CHANNEL2 = CS.FBR.Rx2;
SAT_VEL = CS.GEO.V.V(:);
ALT_RATE = CS.GEO.H_rate(:);
ALTITUDE = CS.GEO.H(:);
   
[SAMPLENO,PULSEPERBURST,BURSTNO,RECNO] = size(CHANNEL1);

TOT_RECORD = RECNO-1;
STR_RECORD = 1;

%DECLARE CONSTANTS:
BLOCK = 20;
ECHO_IN_BURST = 64;
SAMPLE = 4*128;
TRK_WIN = 4*60; %[m]
N = PULSEPERBURST;
%DEFINE: hamming window
n = [1:N]-1;
hamWin = 0.54 - 0.46*cos(2*pi*n/(N-1));
PRF = 18181.82;
%BRI = 11.7e-3; %for SAR mode
BRI = 4*11.7e-3; %for SARIN mode
c = 299792458; 
lam = c/13.575e9; 
k0 = 2*pi/lam;
R = 6378.1370e3; %radius of Earth in meters

%now some rudimentary error checking to make sure your structures are at least
%in the correct format - can't be sure about the data ;)
if (BURSTNO ~= BLOCK)
    BURSTNO
    return;
end
if (PULSEPERBURST ~= ECHO_IN_BURST)
    PULSEPERBURST
    return;
end
if (SAMPLENO ~= SAMPLE)
    SAMPLENO
    return;
end

N = SAMPLENO;
for recN=STR_RECORD:TOT_RECORD
    
    for j=1:BLOCK-1
        for m=1:ECHO_IN_BURST
            chan1 = double(squeeze(CHANNEL1(:,m,j,recN)));
            chan2 = double(squeeze(CHANNEL2(:,m,j,recN)));
            
            chan1Pad = [chan1(1:N/2).' zeros(1,N) chan1(N/2+1:end).'];
            chan1Pad(N/2+1) = chan1(N/2+1)/2;
            chan1Pad(N+N/2+1) = chan1(N/2+1)/2;
            fftChan1 = fft(chan1Pad,N*2)/sqrt(N*2);
            fftChan1 = fftshift(fftChan1);            
            
            chan2Pad = [chan2(1:N/2).' zeros(1,N) chan2(N/2+1:end).'];
            chan2Pad(N/2+1) = chan2(N/2+1)/2;
            chan2Pad(N+N/2+1) = chan2(N/2+1)/2;
            fftChan2 = fft(chan2Pad,N*2)/sqrt(N*2);
            fftChan2 = fftshift(fftChan2);            
            
            %figure;hold on;...
            %     plot(abs(fftChan1),'r');plot(abs(fftChan2),'b');
            
            %form the sum of the power from both channels:
            pow2Chan(m,:) = (abs(fftChan1).^2+abs(fftChan2).^2)/2;
%             figure;hold on;plot(pow2Chan(m,:));
        end
        %the onboard gate retracker is stable within a burst, so can safely
        %add the waveforms within a burst together to get one PLRM (@20Hz):
        PLRM(recN,j,:) = mean(pow2Chan,1); %TODO: should I normalise by the mean ave waveform before finding average?
        %figure;hold on;plot(squeeze(PLRM(recN,j,:)));
    end
    
end            