# nnprime
Simulates precession of neutron mirror neutron system according to the equations derived in E.D. Davis's manuscript. From the theory of mirror neutrons invented by Zurab Berezhiani.

It simulates in either a cylinder, rectangle or sphere. The n-n' mixing phase IS reset on wall collision (alhtough it does this in the code, in reality its not clear to me that this is true, afterall neutrons in interferometry don't consider a wall collision a measurement.)

UCN have a hard v^2 distribution and gravity is accounted for in the starting position via equaiton 4.7 in Bob G's Book. The Fermi potential in the trap is set by the maximum velocity allowed in the trap, there is no consideration for marginal trapping.  

Requires boost libraries. I used 1.64

Makefile requires an export of $BOOST_INCLUDE which is the location of the boost library directory.
I know my makefile is lame, also I know I don't care. If you care, then fix it.

once the thing is made it can be run with: $./runMirror parameterfile outputdata

parameterfile is the name of the parameter file. outputdata is the output data file name

The MirrorParameter.txt file is a good example of a parameter file, also it basically explains everything.

I included a python extraction/plotting script, matlab extraction/plotting is in the ./data/ folder. I will use matlab so it will be more mature in the long run. 

