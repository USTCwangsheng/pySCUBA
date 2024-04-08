## pySCUBA SCUBASketch

1.You can generate structure sketches directly using `gen_alphabetaN_crd.py` and `gen_beta_propellerN_crd.py`. Notice change the program installation path in the files.

2.You can use `SCUBAsketch.py` to generate sketch structures from other coordinates files.

### Tutorial
program name: SCUBAsketch

Function: Generate main chain structure from given secondary structure information. You need to provide: secondary structure type, length range, build direction (n-terminal or C-terminal), starting position.

Run command:
```sh
installation_path/pySCUBA/cpp_bin/SCUBASketch sketch.par 
```

Parameter file:`sketch.par`
```
START SketchPar
SketchFile = 3pvh.txt    #Secondary structure information file.
LinkLoop = 1             #Whether to use Loop to connect the secondary structure. 1 indicates yes and 0 indicates no.
RandomSeed = 123         #RandomSeed.
OutputFile = 3pvh_build  #The prefix of the output file name.
GenerationNumber = 3     #The number of structures generated.
END SketchPar
```

Example of secondary structure information file:`3pvh.txt`
```
H 10 10 N; 10 0 0 ; 0 0 1
E 4 4 C  ; 15 10 0; 0 0 -1
H 20 20 N; 16 20 0; 0 0 1
E 7 7 C  ; 10 10 0; 0 0 -1
H 14 14 N; 0 0 0  ; 0 0 1
E 7 7 C  ; 5 10 0 ; 0 0 -1
E 7 7 N  ; 0 10 0 ; 0 0 1
H 20 20 C; 0 20 0 ; 0 0 -1
```
Each line represents a secondary structure of information.
The first semicolon indicates the type of secondary structure (E for β-strand, H for α-helix), the length range (two integers), the direction of construction (N or C, N for starting from the N terminal, C for starting from the C terminal), and the second semicolon indicates the starting position of the first amino acid Cα coordinate; The second semicolon is followed by the termination direction (vector) of the secondary structure.
The program generates these secondary structures in order.