# plotdraw
A package to help with drawing to plotters.

Basic drawing functions; line, paths, bezier paths, polygon, circle, ellipse

Plotting optimisation; join separate paths with same start/end together to create contiguous paths, convert straight line paths to lines, optimise order of paths for efficient plotting by minimising pen up time

Plotting; output .HPGL and .SVG files

Requires svgwrite for SVG preview functionality.

Check main.py for examples.
