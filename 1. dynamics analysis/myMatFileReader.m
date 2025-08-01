function [GM, Re, degree, C, S] = myMatFileReader(datafile_name)
    % Carica tutte le variabili dal file .mat
    loaded_data = load(datafile_name);

    % Assicurati che le variabili necessarie siano presenti
    required_vars = {'Re', 'GM', 'degree', 'C', 'S'};
    for i = 1:length(required_vars)
        if ~isfield(loaded_data, required_vars{i})
            error('myMatFileReader:MissingVariable', ...
                  'Il file dati non contiene la variabile richiesta: %s', required_vars{i});
        end
    end

    % Estrai i dati
    C = loaded_data.C;
    S = loaded_data.S;
    Re = loaded_data.Re;
    GM = loaded_data.GM;
    degree = loaded_data.degree;
end
