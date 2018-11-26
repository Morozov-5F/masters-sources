settings.prc = false;
size(5cm, 5cm);

import three;

// Set the projection
//currentprojection = orthographic(-3.78563396854908,11.2139874694386,2.4453953248835);
// currentprojection=orthographic(
// camera=(-3.78563396854908,11.2139874694386,2.4453953248835),
// up=(0.000568093225207644,-0.00581802430077633,0.027559505239397),
// target=(-1.77635683940025e-15,-8.88178419700125e-16,0),
// zoom=0.47402770397409);

currentprojection=orthographic(
camera=(-6.51556088388528,10.6713482986402,2.40117549478684),
up=(0.00181379536249954,-0.00585051069210866,0.0309226533940109),
target=(0,-1.77635683940025e-15,-4.44089209850063e-16),
zoom=0.47402770397409);

// Draw coordinates
// draw(O -- 2X, L = Label("$x$", position = EndPoint));
// draw(O -- 2Y, L = Label("$y$", position = EndPoint));
// draw(O -- 2Z, L = Label("$z$", position = EndPoint));

// Draw foundation
triple M0 = (2.5, -1.5, 2);
triple M1 = (0, 0, 0);
triple M2 = (4 ,0 ,0);
triple M3 = (3, -3, 0);
triple M4 = (2.25, 0, 0);

dot(M0, L = Label("$M_0$", align = Relative(N)));
dot(M1, L = Label("$M_1$", align = Relative(S)));
dot(M2, L = Label("$M_2$", align = LeftSide));
dot(M3, L = Label("$M_3$", align = Relative(S), position = MidPoint * 0.8));
dot(M4, L = Label("$M_4$", align = Relative(S)));

draw(M1 -- M2);
draw(M1 -- M3, dashed);
draw(M2 -- M3, dashed);

draw(M0 -- M1);
draw(M0 -- M2);
draw(M0 -- M3, dashed);
draw(M0 -- M4);

guide3 M04 = M4 -- M0;
guide3 M41 = M4 -- M1;
draw(arc(M4, point(M04, 0.1), point(M41, 0.1)), L = scale(0.5) * Label("$\varphi$", position = MidPoint, align = Relative(E)));


guide3 M01 = M0 -- M1;
guide3 M02 = M0 -- M2;
guide3 M04 = M0 -- M4;

draw(arc(M0, point(M01, 0.3), point(M04, 0.3)), L = scale(0.5) * Label("$\alpha_{14}$", align = Relative(S)));
draw(arc(M0, point(M02, 0.28), point(M04, 0.28)), L = scale(0.5) * Label("$\alpha_{24}$", align = Relative(S)));

