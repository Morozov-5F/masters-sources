settings.prc = false;
size(5cm, 5cm);

import three;

// Set the projection
currentprojection = orthographic(-10, 6, 4);

// Draw coordinates
// draw(O -- 2X, L = Label("$x$", position = EndPoint));
// draw(O -- 2Y, L = Label("$y$", position = EndPoint));
// draw(O -- 2Z, L = Label("$z$", position = EndPoint));

// Draw foundation
triple M0 = (2.5, -1.5, 2);
triple M1 = (0, 0, 0);
triple M2 = (4 ,0 ,0);
triple M3 = (3, -3, 0);

dot(M0, L = Label("$M_0$", align = Relative(N)));
dot(M1, L = Label("$M_1$", align = Relative(S)));
dot(M2, L = Label("$M_2$", align = LeftSide));
dot(M3, L = Label("$M_3$"));

draw(M1 -- M2, L = scale(0.5) * Label("$d_{12}$", align = LeftSide));
draw(M1 -- M3, L = scale(0.5) * Label("$d_{13}$"));
draw(M2 -- M3, dashed, L = scale(0.5) * Label("$d_{23}$", align = Relative(S), position = MidPoint));

draw(M0 -- M1, L = scale(0.5) * Label("$l_1$", align = RightSide));
draw(M0 -- M2, L = scale(0.5) * Label("$l_2$", align = LeftSide));
draw(M0 -- M3, L = scale(0.5) * Label("$l_3$", align = RightSide));

guide3 M01 = M0 -- M1;
guide3 M02 = M0 -- M2;

draw(arc(M0, point(M01, 0.1), point(M02, 0.1)), L = scale(0.5) * Label("$\alpha_{12}$", align = Relative(S)));
