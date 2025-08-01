function plot_macrofeature_scores(score_struct)
% Raggruppa e plottia i punteggi ottenuti per ciascuna macro-feature (cioè id comune)

% Estrai gli indici delle facce con punteggio
valid_idx = find([score_struct.score] > 0);

% Estrai tutti gli ID unici (macro-feature)
all_ids = arrayfun(@(s) s.id, score_struct(valid_idx));
unique_ids = unique(all_ids);

% Inizializza vettori
num_features = numel(unique_ids);
obtained_scores = zeros(1, num_features);
total_scores = zeros(1, num_features);

% Per ogni macro-feature calcola il punteggio ottenuto e il totale disponibile
for i = 1:num_features
    current_id = unique_ids(i);
    
    idx = valid_idx(all_ids == current_id);
    
    for j = idx
        completed = score_struct(j).completeness;
        s = score_struct(j).score;
        obtained_scores(i) = obtained_scores(i) + s * sum(completed) / length(completed);
        total_scores(i) = total_scores(i) + s;
    end
end

% Percentuale ottenuta
percentage = 100 * obtained_scores ./ total_scores;

% Plot
figure;
bar(percentage);
xlabel('Macro-feature ID');
ylabel('Completion [%]');
title('Score completion per macro-feature');
xticks(1:num_features);
xticklabels(string(unique_ids));
grid on;
ylim([0, 110]);

end
