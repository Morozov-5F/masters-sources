//import three;
import solids;
import three;

unitsize(2.5inch);
//size(5inch);

settings.render=0;
settings.prc=false;

currentprojection=orthographic(3, 1, 1);
currentlight=nolight;

real RE=0.6, RS=0.7, inc=100, lat=45, lon=45, tlat=50, tlon=100;

//Added specification for mesh pen
// draw(surface(Earth), surfacepen=lightgrey+opacity(.6), meshpen=0.6*white);
//The next line is no longer necessary.
//draw(Earth,m=10,gray);

real R = 0.7, r = 1;
triple center =  (R, 0, 0);
int nu = 36, nv = 36;

path3 crossSection = Circle(r = r, c = center, normal = Y, n = nu);

surface torus = surface(crossSection, c=(0,0,0), axis=Z, n=nv, angle1=0, angle2=360);

draw(torus, surfacepen=grey+opacity(.6), meshpen=0.4*black+opacity(0.7));

triple Mi = (0, 0,  0.714143);
triple Mj = (0, 0, -0.714143);

real theta = 30;
real phi = 60;

triple M0 = ((R + r * Cos(theta)) * Cos(phi),
             (R + r * Cos(theta)) * Sin(phi),
              r * Sin(theta));

draw(Mi -- Mj, p=linewidth(1), L = Label("$d_{12}$", align = Relative(SW), p = fontsize(20pt)));
draw(Mi -- M0, p=linewidth(1));
draw(Mj -- M0, p=linewidth(1));

dot(Mi, L = Label("$M_1$", align = Relative(SW), p = fontsize(20pt)));
dot(Mj, L = Label("$M_2$", align = Relative(SW), p = fontsize(20pt)));
dot(M0, L = Label("$M_0$", align = Relative(SE), p = fontsize(20pt)));