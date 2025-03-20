# Processing Bezier curves for machine learning
The main file is
```
processCurves.ipynb
```
It contains 3 main algorithms to process the curves :
1. Uniformal law
2. Geometrical law
3. Curvature based

## Known issues
- The curvature is too small for lines, sometimes 0, so delta always returns 0. The problem was fixed by separating the points, but that's a temporary fix.
- The documentation is lacking 😕
- Too much DRY and unoptimised code. But the main goal is for the program to work for all curvatures first and then optimize.

## Work to do
- Test the curvature based algorithm on bezier curves only (a single one, then multiple).
- Find a work-around for lines and zero/infinite derivatives in the intersection of lines.
- Add optimised laws (maybe stochastic ones)
- Test for all curves provided by Louis-Alexis.