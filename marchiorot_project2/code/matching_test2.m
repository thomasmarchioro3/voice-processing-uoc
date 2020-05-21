rng(1184445);

clc, close, clear;
Ns = 1000; % number of samples
Fs = 16000;
Delta = 40; % threshold for matching

n1 = 20; % number of peaks for the beginning of the frame
n2 = 20; % number of peaks for the end of the frame
A_max = 100;

F_1 = rand(1, n1)*Fs; F_1 = sort(F_1);
F_2 = rand(1, n2)*Fs; F_2  = sort(F_2);
AMP_1 = randi(A_max, 1, n1);
AMP_2 = randi(A_max, 1, n2);
PH_1 = 2*pi*rand(1, n1);
PH_2 = 2*pi*rand(1, n2);

L_1 = length(F_1);
L_2 = length(F_2);

[idx_pairs1, L_new1] = findMatchesM(F_1, F_2, L_1, L_2, Delta);
[idx_pairs2, L_new2] = findMatchesT(F_1, F_2, L_1, L_2, Delta);




return;

AMP_1_new = zeros(L_new,1);
AMP_2_new = zeros(L_new,1);
PH_1_new = zeros(L_new,1);
PH_2_new = zeros(L_new,1);
F_1_new = zeros(L_new,1);
F_2_new = zeros(L_new,1);

for k = 1:L_new
    if idx_pairs(1, k) == 0
        % birth
        %plot([1, 2], [F_2(idx_pairs(2, k)), F_2(idx_pairs(2, k))], 'go-');
        F_1_new(k) = F_2(idx_pairs(2, k));
        F_2_new(k) = F_2(idx_pairs(2, k));
        AMP_1_new(k) = AMP_2(idx_pairs(2, k));
        AMP_2_new(k) = AMP_2(idx_pairs(2, k));
        PH_1_new(k) = PH_2(idx_pairs(2, k));
        PH_2_new(k) = PH_2(idx_pairs(2, k));
    elseif idx_pairs(2, k) == 0
        % death
        %plot([1, 2], [F_1(idx_pairs(1, k)), F_1(idx_pairs(1, k))], 'ko-');
        F_1_new(k) = F_1(idx_pairs(1, k));
        F_2_new(k) = F_1(idx_pairs(1, k));
        AMP_1_new(k) = AMP_1(idx_pairs(1, k));
        AMP_2_new(k) = AMP_1(idx_pairs(1, k));
        PH_1_new(k) = PH_1(idx_pairs(1, k));
        PH_2_new(k) = PH_1(idx_pairs(1, k));
    else 
        % match
        %plot([1, 2], [F_1(idx_pairs(1, k)), F_2(idx_pairs(2, k))], 'ro-');
        F_1_new(k) = F_1(idx_pairs(1, k));
        F_2_new(k) = F_2(idx_pairs(2, k));
        AMP_1_new(k) = AMP_1(idx_pairs(1, k));
        AMP_2_new(k) = AMP_2(idx_pairs(2, k));
        PH_1_new(k) = PH_1(idx_pairs(1, k));
        PH_2_new(k) = PH_2(idx_pairs(2, k));
    end
end

figure(1);
for k = 1:L_new
    %disp(idx_pairs(:, k));
    if idx_pairs(1, k) == 0
        % birth
        plot([1, 2], [F_2(idx_pairs(2, k)), F_2(idx_pairs(2, k))], 'go-');
    elseif idx_pairs(2, k) == 0
        % death
        plot([1, 2], [F_1(idx_pairs(1, k)), F_1(idx_pairs(1, k))], 'ko-');
    else 
        % match
        plot([1, 2], [F_1(idx_pairs(1, k)), F_2(idx_pairs(2, k))], 'ro-');
    end
    hold on;
end
hold off;
xlim([0.5, 2.5]);

figure(2);
for k = 1:L_new
    plot([1, 2], [F_1_new(k), F_2_new(k)], 'ko-');
    hold on;
end
hold off;
xlim([0.5, 2.5]);


function [idx_pairs, L_new] = findMatchesT(F_1, F_2, L_1, L_2, Delta)

  L_new = 0;

  idx_pairs = zeros(2, L_1 + L_2);

  i = 1; j = 1;
  while i <= L_1 && j <= L_2
      if F_2(j) < F_1(i) - Delta
          % birth
          L_new = L_new + 1;
          idx_pairs(1, L_new) = 0;
          idx_pairs(2, L_new) = j;
          j = j + 1;
      elseif F_2(j) > F_1(i) + Delta
          % death
          L_new = L_new + 1;
          idx_pairs(1, L_new) = i;
          idx_pairs(2, L_new) = 0;
          i = i + 1;
      else
        % match
        L_new = L_new + 1;
        cond = true;
        while cond && j < L_2
            cond = false;
            improvement = abs(F_1(i) - F_2(j)) > abs(F_1(i) - F_2(j + 1));
            if improvement
                % birth
                idx_pairs(1, L_new) = 0;
                idx_pairs(2, L_new) = j;
                j = j + 1;
                L_new = L_new + 1;
                % continue cycling
                cond = true;
            end
        end
        idx_pairs(1, L_new) = i;
        idx_pairs(2, L_new) = j;
        i = i + 1; j = j + 1;
      end
        
  end

  idx_pairs = idx_pairs(:, 1:L_new);

end

function [idx_pairs, L_new] = findMatchesM(F_1, F_2, L_1, L_2, Delta)
  IndexMatch = zeros(L_1 + L_2, 2);
  i = 1; j = 1; k = 1;
  while (i <= L_1) && (j <= L_2)

      if F_2(j) > F_1(i) + Delta %THANATOS

          IndexMatch(k,1) = i;
          IndexMatch(k,2) = 0;

          i = i + 1;

      elseif F_2(j) < F_1(i) - Delta %GENESIS

          IndexMatch(k,1) = 0;
          IndexMatch(k,2) = j;

          j = j + 1;

      else

          if IndexMatch(k,2) == 0 %NO CANDIDATE STORED YET

              IndexMatch(k,1) = i;
              IndexMatch(k,2) = j;

          elseif abs(F_1(i) - F_2(j-1)) < abs(F_1(i) - F_2(j)) %WORSE CANDIDATE => CONTINUE

               i = i + 1;
               k = k - 1;

          else %BETTER CANDITATE => + GENESIS TO PREVIOUS CANDIDATE

              IndexMatch(k-1,1) = 0;

              IndexMatch(k,1) = i;
              IndexMatch(k,2) = j;

          end

          j = j + 1;

      end

      k = k + 1;

  end

  k = k - 1;
  L_new = k;

  IndexMatch = IndexMatch(1:1:k, 1:1:2);

  idx_pairs = IndexMatch';

end




