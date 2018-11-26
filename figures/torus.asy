//import three;
import solids;
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

real r1 = 0.7, r0 = 1;
int nu = 36, nv = 36;

path3 crossSection = Circle(r=r0, c=(r1,0,0), normal=Y, n= nu);

surface torus = surface(crossSection, c=(0,0,0), axis=Z, n=nv, angle1=0, angle2=360);

draw(torus, surfacepen=grey+opacity(.5), meshpen=0.6*black+opacity(0.9));

triple Mi = (0, 0, 0.714143);
triple Mj = (0, 0, -0.714143);

draw(Mi -- Mj, p=linewidth(1));

dot(Mi);
dot(Mj);