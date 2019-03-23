    # function [cD cA] = getDWT(X,N,Name)
    
    # finds the discrete wavelet transform at level N for signal X using the
# wavelet specified by Name.
    
    ## Inputs:
# X: the original signal
# N: the decomposition level
# Name: the wavelet name to use
    
    ## Outputs:
# cD is a N-row matrix containing the detail coefficients up to N levels
# cA is the same for the approximations
    
    # This code was developed by David Springer for comparison purposes in the
# paper:
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
def getDWT(X=None,N=None,Name=None,*args,**kwargs):
    varargin = getDWT.varargin
    nargin = getDWT.nargin

    #No DWT available for Morlet - therefore perform CWT:
    if (strcmp(Name,'morl')):
        c=cwt(X,arange(1,N),'morl')

        cD=copy(c)

        cA=copy(c)

    else:
        #Preform wavelet decomposition
        c,l=wavedec(X,N,Name,nargout=2)

        #decomposition 
        len_=length(X)

        cD=zeros(N,len_)

        for k in arange(1,N).reshape(-1):
            d=detcoef(c,l,k)

            d=ravel(d).T

            d=d(ones(1,2 ** k),arange())

            cD[k,arange()]=wkeep1(ravel(d).T,len_)

        cD=ravel(cD)

        I=find(abs(cD) < sqrt(eps))

        cD[I]=zeros(size(I))

        cD=reshape(cD,N,len_)

        #Reorder the approximations based on the structure of the wavelet
    #decomposition 
        len_=length(X)

        cA=zeros(N,len_)

        for k in arange(1,N).reshape(-1):
            a=appcoef(c,l,Name,k)

            a=ravel(a).T

            a=a(ones(1,2 ** k),arange())

            cA[k,arange()]=wkeep1(ravel(a).T,len_)

        cA=ravel(cA)

        I=find(abs(cA) < sqrt(eps))

        cA[I]=zeros(size(I))

        cA=reshape(cA,N,len_)

    