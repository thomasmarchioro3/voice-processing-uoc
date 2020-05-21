function get_data(folder, N_mel, Horizon, Shift)
if nargin == 0
    folder = 'test';
end
if nargin <= 1
    N_mel = 25;
end
if nargin <= 2
    Horizon = 20;
end
if nargin <= 3
    Shift = 5;
end

addpath(genpath('sap-voicebox'));
path = strcat('TIMIT/', folder);
addpath(genpath(path));

dirs = dir(path);
dirs = dirs(3:end); % remove '.' and '..'
data = [];
for ii = 1:length(dirs)
    currdir = dirs(ii);
    speaker.name = currdir.name;
    files = dir(strcat(currdir.folder, '/', currdir.name));
    files = files(3:end); % remove '.' and '..'    
    melcoeff = [];
    for jj = 1:length(files)
        file = files(jj);
        file = strcat(file.folder, '/', file.name);
        [sig, Fs] = audioread(file);
        Win = 'N';
        p = floor(3*log(Fs));
        Hor = Horizon*Fs/1000;
        Sh = Shift*Fs/1000;
        % compute the Mel Cepstrum Coefficients for jj-th file
        coeff = v_melcepst(sig, Fs, Win, N_mel, p, Hor, Sh);
        % append the coefficients
        melcoeff = [melcoeff; coeff]; 
    end    
    speaker.coeff = melcoeff;
    data = [data, speaker];
end

outfile = strcat(folder, '_data.mat');
save(outfile, 'data');
clear all;
end