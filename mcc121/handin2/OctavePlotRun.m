[hFig, hSurf] = OctavePlot(load('X.txt'), "X", load('Y.txt'), "Y", load('Z.txt'), "Z", 1);
 while(ishandle(hFig)) # sleep until figure is closed
 sleep(0.5);
 end
 exit(0) # close octave and free up some memory

