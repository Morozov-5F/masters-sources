settings.outformat="png";
settings.prc = false;

import three;

size(4cm, 0);

draw(-2X -- 2X, blue, arrow=Arrow3(emissive(blue)), L=Label("$x'$", position=EndPoint, align=W, black));  //x-axis
draw(-2Y -- 2Y, green, arrow=Arrow3(emissive(green)), L=Label("$y'$", position=EndPoint, black)); //y-axis
draw(-2Z -- 1Z, red, arrow=Arrow3(emissive(red)), L=Label("$z'$", position=EndPoint, black));   //z-axis

triple M = (1, 1.5, 0);
triple PO = M - 1Z;

draw(PO--(0, 0, 0), black, arrow=Arrow3(HookHead2, emissive(black), size=0.9));
draw((0, 0, 0)--M, dashed);
draw(PO -- M, dashed);

dot(PO, L = Label("$PO$"));
dot(M, L = Label("$M$"));

draw(arc((0, 0, 0), point(X, 0.1), point(M, 0.1)), L = scale(0.8) * Label("$\alpha$", align = Relative(S)), arrow= Arrow3(HookHead2, emissive(black), size=1.0));

draw(arc((0, 0, 0), point(M * 0.5, 0.1), point(PO * 0.5, 0.1)), L = scale(0.7) * Label("$\varepsilon$", align = Relative(SE)), arrow= Arrow3(HookHead2, emissive(black), size=1.0));
