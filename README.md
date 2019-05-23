# plotdraw
A package to help with drawing to plotters.

Basic drawing functions; line, paths, bezier, cubic bezier, polygon, rounded polygon, circle, ellipse, arcs

Plotting optimisation; automatically join all separate paths with same start/end together to create contiguous paths, convert straight line paths to lines, optimise order of paths for efficient plotting by minimising pen up time

Plotting; output .HPGL and .SVG files

Requires svgwrite for SVG preview functionality.

Check main.py for examples.
