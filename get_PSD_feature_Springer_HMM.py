    #cfunction [psd] = get_PSD_feature_Springer_HMM(data, sampling_frequency, frequency_limit_low, frequency_limit_high, figures)
    
    # PSD-based feature extraction for heart sound segmentation.
    
    ## INPUTS:
# data: this is the audio waveform
# sampling_frequency is self-explanatory
# frequency_limit_low is the lower-bound on the frequency range you want to
# analyse
# frequency_limit_high is the upper-bound on the frequency range
# figures: (optional) boolean variable to display figures
    
    ## OUTPUTS:
# psd is the array of maximum PSD values between the max and min limits,
# resampled to the same size as the original data.
    
    # This code was developed by David Springer in the paper:
# D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
# Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
    
    ## Copyright (C) 2016  David Springer
# dave.springer@gmail.com
    
    # This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
    
    # This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
    
    # You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    
@function
def get_PSD_feature_Springer_HMM(data=None,sampling_frequency=None,frequency_limit_low=None,frequency_limit_high=None,figures=None,*args,**kwargs):
    varargin = get_PSD_feature_Springer_HMM.varargin
    nargin = get_PSD_feature_Springer_HMM.nargin

    if nargin < 5:
        figures=0

    
    # Find the spectrogram of the signal:
    __,F,T,P=spectrogram(data,sampling_frequency / 40,round(sampling_frequency / 80),arange(1,round(sampling_frequency / 2),1),sampling_frequency,nargout=4)

    if (figures):
        figure()
        surf(T,F,dot(10,log(P)),'edgecolor','none')
        axis('tight')
        view(0,90)
        xlabel('Time (Seconds)')
        ylabel('Hz')
        pause()
    
    __,low_limit_position=min(abs(F - frequency_limit_low),nargout=2)

    __,high_limit_position=min(abs(F - frequency_limit_high),nargout=2)

    # Find the mean PSD over the frequency range of interest:
    psd=mean(P(arange(low_limit_position,high_limit_position),arange()))

    if (figures):
        t4=(arange(1,length(psd))) / sampling_frequency

        t3=(arange(1,length(data))) / sampling_frequency

        figure('Name','PSD Feature')
        plot(t3,(data - mean(data)) / std(data),'c')
        hold('on')
        plot(t4,(psd - mean(psd)) / std(psd),'k')
        pause()
    