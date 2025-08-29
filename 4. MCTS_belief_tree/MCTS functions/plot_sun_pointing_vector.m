function plot_sun_pointing_vector(t_model, frame)

cmap = [linspace(1, 1, length(t_model))', linspace(1, 0.7, length(t_model))', linspace(0.5, 0, length(t_model))'];
for i = 1:length(t_model)
    r_sun = cspice_spkgeo(2000433, t_model(i), 'ECLIPJ2000', 10);
    r_sun = r_sun(1:3);
    eclip2eros = cspice_pxform('ECLIPJ2000', frame, t_model(i));
    r_sun = eclip2eros*(-r_sun);
    r_sun = r_sun/norm(r_sun)*40;
    plot3([0, r_sun(1)], [0, r_sun(2)], [0, r_sun(3)], 'Color', cmap(i,:));
end
