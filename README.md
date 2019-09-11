# nnprime
Simulates precession of neutron mirror neutron system according to the equations derived in E.D. Davis's manuscript. From the theory of mirror neutrons invented by Zurab Berezhiani.

It simulates in either a cylinder, rectangle or sphere. The n-n' mixing phase IS reset on wall collision (alhtough its not clear to me this is true, afterall neutrons in interferometry don't consider a wall collision a measurement.)

UCN have a hard v^2 distribution and gravity is accounted for in the starting position via equaiton 4.7 in Bob G's Book. The Fermi potential in the trap is set by the maximum velocity allowed in the trap, there is no consideration for marginal trapping.  

Requires boost libraries. I used 1.64

Makefile requires an export of $BOOST_INCLUDE which is the location of the boost library directory.
I know my makefile is lame, also I know I don't care. If you care, then fix it.

once the thing is made it can be run with. 

./runMirror parameterfile outputdata

paremeterfile is the name of the parameter file, if its not obvious. 
outputdata file is the outputdatafile, i really hope you figure this out on your own. 

The MirrorParameter.txt file is a good example of a parameter file, also it basically explains everything.

also I think there needs to be a ./data/ that git is ignoring. that is where the data is saved.

I included a python script to plot data, matlab script will come soon. that is what I will use so it will be more mature in the long run. 

