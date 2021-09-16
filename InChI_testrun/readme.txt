Test runs on DataWarrior's libraries of random molecules (n = 1k, n = 5k)
written into a .sdf file and subsequently assigned by inchi-1 their
InChI and InChiKey.

inchi-1 is version 1.06 of the InChI Trust's reference implementation,
compiled for Linux, version 1.06 (December 2020).  Running md5sum yields

d460208e006c244778e8f772d6f99d66  inchi-1

After provision of the executable bit, the program was run in the pattern
of

./inchi-1 test.sdf -AuxNone -Tabbed -Key

to read the content of file test.sdf and to assign the InChI.  The used
options instruct the program further to

+ skip the auxiliary layer in the InChi (-AuxNone)
+ report the results as a tabulator separated table (the compound counter
  hereby is the first column, followed by the InChI)
+ to append per entry the InChIKey (content of the third column).


Files created retain the root of the name of the input file.  By file
extension

+ .log retains the comments sent to the CLI
+ .prb possibly would report if there is some significant problem for
  one entry processed
+ .txt retains the permanent record of the newly assigned data.

END
