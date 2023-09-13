% Known time points
known_times = [0, 1/3, 2/3, 1];
	
	
	


% Known x and y coordinates, replace x0, y0, ..., x3, y3 with real values
known_xs = [0, -.7, -0.05, -0.65];
known_ys = [0, -.35, -.53, -0.7];

% Times at which to interpolate
inter_times = linspace(0, 1, 1024);

% Perform the interpolation
inter_xs = interp1(known_times, known_xs, inter_times);
inter_ys = interp1(known_times, known_ys, inter_times);

% Combine the interpolated x and y coordinates
inter_coords = [inter_xs; inter_ys].';

error_fusion = sqrt(sum((inter_coords - array2).^2, 2));