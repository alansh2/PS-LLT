function initutils
%INITUTILS
filename = mfilename('fullpath');
filepath = fileparts( filename );

addpath([filepath,'/core']);
addpath([filepath,'/custom']);
addpath([filepath,'/poly-test']);