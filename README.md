# sc_golf
This is a prototype for a game based on straightedge and compass constructions.  
(The objective is to approximate a given set of points as closely as possible starting from another set, using a fixed set of moves.)

REQUIREMENTS:

- basics.c/h (from whw_clib)
- system.c/h (from whw_clib)
- SDL2

To define a circle:
  - While pressing and holding 'c', click near a pre-defined circular vertex (not the purple diamond ones) to select the center.
  - Release 'c', and click near another circular vertex to select a radial point.
To define a line:
  - Essentially the same process as drawing a circle except with 'l' instead of 'c' (two points define a line.)

To define a point:
- Press 'p' (and release) to enter 'point-adding mode' (you can exit this at any time by pressing 'esc'.)
- To select a curve, press 'c' or 'l' to enter 'curve-selecting mode' ('c' for circles and 'l' for lines.)
- Use the 'up' and 'down' arrows to move between different curves, and press 'enter' once you have chosen a curve.
- After selecting the first curve, press 'c' or 'l' again to select the second, and repeat the above process.
- If you chose curves that intersect at two points, use the up and down keys to select an intersection point.
- Press 'enter' until the new point has been added and the 'highlights' disappear (this part definitely needs to be improved.)
- You can also try clicking near a line or circle while in selection mode (but I haven't tested/ironed it out yet.)

To undo:
- Press (and release) 'ctrl', and then press 'u' however many times you need (this needs to be adjusted.)

Scoring: 
I was a bit torn between a basic least-squares score (which penalizes 'favoritism' and inaccuracy and rewards equanimity) and a 
fractional lp score (which penalizes inexactness/imprecision.) Currently I'm experimenting with a l(1/4) score, but that could 
easily change in the next version (and of course you are welcome to adjust it to suit your own preferences.) 
I also thought it might be interesting to keep individuated tallies of each elementary operation (including 'undos' and
adding points), and to incorporate the tally in the score. Another version might instead have pre-defined neighborhoods for 
each of the 'holes' that need to be reached.

Update 11.4.2023
In the most recent uploaded version of sc_golf.c, the score is based on an "electrostatic" potential, where target points
are negatively charged and the constructed points are positive. The game registers if all points have been reached to within
a prescribed tolerance, but I have yet to implement a way to adjust the tolerance (or "hole width", in golf terms.) 
In the next iteration, I hope to incorporate features to (a) change the game settings and scoring, (b) save and resume progress,
and (c) analyze data (e.g. statistical correlations between (randomly chosen) hole positions and various scores.
