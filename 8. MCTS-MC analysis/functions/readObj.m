function [obj] = readObj(filename)
    fid = fopen(filename);
    obj.v = [];
    obj.f.v = [];

    while ~feof(fid)
        tline = strtrim(fgetl(fid));
        if startsWith(tline, 'v ')
            coords = sscanf(tline(3:end), '%f %f %f');
            obj.v(end+1,:) = coords';
        elseif startsWith(tline, 'f ')
            parts = strsplit(tline(3:end));
            face = zeros(1, length(parts));
            for i = 1:length(parts)
                % Estrarre indice del vertice prima dello slash (formato: v, v/vt, v//vn, v/vt/vn)
                slash_pos = strfind(parts{i}, '/');
                if isempty(slash_pos)
                    face(i) = str2double(parts{i});
                else
                    face(i) = str2double(parts{i}(1:slash_pos(1)-1));
                end
            end
            obj.f.v(end+1,:) = face;
        end
    end
    fclose(fid);
end
