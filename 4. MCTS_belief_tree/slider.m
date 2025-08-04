function slider(F, V, new_scores, new_known_map, XX, TT)
    N = size(XX,1);

    figure('Name', 'Animazione con slider', 'Position', [100 100 800 600]);
    ax = axes('Position', [0.1 0.2 0.8 0.75]);
    hold(ax, 'on'); axis(ax, 'equal'); view(3);

    % Plotta l’asteroide (mesh)
    plotEllipsoidWithKnownRegion(F, V, new_scores, new_known_map);

    % Orbita completa in azzurro chiaro
    plot3(ax, XX(:,1), XX(:,2), XX(:,3), 'Color', [0.6 0.8 1], 'LineWidth', 1.5);

    % Handle per probe
    h_probe = plot3(ax, XX(1,1), XX(1,2), XX(1,3), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

    % Inizializza handle del vettore Sole
    h_sun = plot_sun_pointing_vector(TT(1)); % <-- Modifica la funzione per restituire l'handle!

    % Slider
    slider = uicontrol('Style','slider', ...
        'Min',1, 'Max',N, 'Value',1, ...
        'SliderStep', [1/(N-1)  10/(N-1)], ...
        'Position',[120 30 560 25], ...
        'Callback',@slider_callback);

    % Label frame
    frame_label = uicontrol('Style','text','Position',[360 5 80 20],...
        'String','Frame: 1','FontSize',12);

    function slider_callback(~,~)
        idx = round(get(slider, 'Value'));
        set(h_probe, 'XData', XX(idx,1), 'YData', XX(idx,2), 'ZData', XX(idx,3));
        if isvalid(h_sun), delete(h_sun); end % Cancella il vecchio vettore Sole
        h_sun = plot_sun_pointing_vector(TT(idx)); % Plotta il nuovo vettore Sole
        set(frame_label, 'String', sprintf('Frame: %d', idx));
        drawnow;
    end

function h = plot_sun_pointing_vector(t)
    % Calcola vettore Sole per il tempo t (scalato di lunghezza fissa)
    r_sun = cspice_spkgeo(2000433, t, 'ECLIPJ2000', 10);
    r_sun = r_sun(1:3);
    eclip2eros = cspice_pxform('ECLIPJ2000', 'IAU_EROS', t);
    r_sun = eclip2eros*(-r_sun);
    r_sun = r_sun/norm(r_sun)*40; % adatta la lunghezza del vettore se vuoi

    % Colore giallo/arancio fisso per il Sole
    colore = [1, 0.7, 0];
    % Plotta e restituisci l'handle
    h = plot3([0, r_sun(1)], [0, r_sun(2)], [0, r_sun(3)], 'Color', colore, 'LineWidth', 2);
end

end