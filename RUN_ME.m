function RUN_ME
% Figures for Moran, Tikhonov
% Code by Jacob Moran and Mikhail Tikhonov
%
% General comment:
%   * Code uses "diversity" interchangeably with "richness" in phrasing
%   
% Potential issues and suggested corrections:
%   *   If an error is thrown from function ukmeans.m, it is most likely
%   due to a MATLAB version issue. Just simply replace ukmeans with
%   MATLAB's in-built kmeans function and results will be similar.
%
% Will reuse precomputed data from the data folder if available
% Remove a data file form that folder to recompute it from scratch
%
% Uses:
%   myColorMap - by Dr Ahmed Abass (c) 2016
%   slightly modifies kmeans by The MathWorks, Inc. (c) 1993-2018 to make
%   new function ukmeans (uncertain k-means algorithm with weighted
%   distance measure w.r.t. uncertainty in points to be clustered)
%
% License conditions applicable to the above:
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

num_invaders = 10;
tStart = tic;
if ~exist('figureData','dir')
    mkdir('figureData');
end
tic
makeFigure2(num_invaders);
toc
tic
makeFigure3;
toc
tic
makeFigure4;
toc
tic
makeFigureNclasses;
toc
tic
makeFigureGLVinvaders(num_invaders);
toc
tic
makeFigureCRMinvaders(num_invaders);
toc
tic
makeFigureNonStandardizedInputs;
toc
tic
makeFigureODcontrol;
toc
tic
makeFigureDownSample;
toc
tic
makeFigureRelabelTest;
toc
tic
makeFigureRandomizeX;
toc
tic
makeFigurePermuteY;
toc
fprintf('Done.\n');
toc(tStart)
end