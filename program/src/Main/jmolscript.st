zap;
load "mpbecinput.pqr";
isoSurface yy ignore {solvent} sasurface 1.4;
write isosurface "jmolsurface.jvxl";
zap;
load "mpbecinput.pqr";
isoSurface "jmolsurface.jvxl" COLORSCHEME rwb MAP "MEP.dx" translucent 0.0;
write zip "jmol_visualization_files.zip"
