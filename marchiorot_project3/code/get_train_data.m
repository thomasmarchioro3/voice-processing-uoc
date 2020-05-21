clc, close, clear;

Horizon  = 30;   %30ms - window length
OrderLPC = 6;   %order of LPC
data_dir = 'TIMIT/train';
D = [pwd,'/',data_dir];
Gmax = 0;
g_array = [];
T = dir(fullfile(D, '**/*.*'));  %get list of files and folders in any subfolder
N = {T(~[T.isdir]).folder}; % files in subfolder.
C = {T(~[T.isdir]).name}; % files in subfolder.
for jj = 1:numel(C)
    F = fullfile(N{jj},C{jj});
    if strcmp(F(end-3:end), '.wav') == 1
        [G, gi, ~] = get_LPC_params(F, OrderLPC, Horizon);
        Gmax = max(Gmax, max(G));
        g_array = [g_array, gi];
    end
end

save('train_data.mat', 'Gmax', 'g_array');