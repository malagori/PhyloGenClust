PhyloGenClust (Phylogenetic Clustering of Gene Families)

Installation guide to PhyloGenClust:

	Required packages:
	a) pip install numpy 
	b) pip install DendroPy
	c) pip install biopython
	d) pip install setuptools
	e) Muscle
	f) FastTree 2.1.7 (user is free to use any FastTree version) 
	g) Notung-2.6 (user is free to use any Notung version)
	h) EMBOSS-6.3.1 for PHYLIPNEW-3.69.650 (for ubuntu)
		sudo apt-get install embassy-phylip 
	   	
	   EMBOSS-6.3.1 installation on Mac osx:
        
        Download EMBOSS-6.3.1 from http://emboss.sourceforge.net/download/
        
        1) tar -zxvf EMBOSS-6.x.x.tar.gz 
        2) cd EMBOSS-6.x.x
        3) ./configure --prefix=/bubo/home/h1/mehmood/install-emboss/
        4) make
		
        Download PHYLIPNEW from http://emboss.sourceforge.net/apps/release/6.6/embassy/index.html
		5) mkdir embassy
		6) mkdir install-emboss
		7) cd embassy
		8) tar -zxvf phylipnew.x.x.tar.gz
		9) cd phylipnew
		10) ./configure --prefix=/path-to-folder/install-emboss/
		11) make 
		12) make install
        (for more instructions, see : http://permalink.gmane.org/gmane.science.biology.emboss/636)
	e) RapidNJ (Optional)
	
Export path to external libraries:
		
	Muscle path:
			export PATH=/path-to-muscle-package/:$PATH
			
	FastTree path:
			export PATH=/path-to-fasttree-package/:$PATH
			
	Notung path:
			export PATH=/path-to-Notung-package/:$PATH

Run GFD:

	>runGFD --help
	usage: runGFD [-h] [-d seqFile] [-st sTree] [-f format] [-t seqType]
              [-nomsa nomsa] [-i interleaved] [-b reps] [-z zeta]
              [-m rsMethod] [-s seed] [-l loadSeed] [-o outFile] [-od outDir]
              [-np NP] [-rapidnj RAPIDNJ] [-fastest FASTEST] [-gtr GTR]
              [-wag WAG] [-gamma GAMMA]

	Parse input arguments and print output.

	optional arguments:
 	-h, --help        show this help message and exit

  	-d seqFile        Specify path to the sequence file

  	-st sTree         Specify path to the species tree file

 	-f format         Specify format of sequence file: F (Fasta); P (Phylip); G
                    (Genbank)

  	-t seqType        Specify type of data: d (dna); p (protein); r (rna)
                    (default='d')

  	-nomsa nomsa      Generate MSA. Use this option, If the input sequences are
                    not aligned

 	-i interleaved    Specify if the sequence data is interleaved
                    (default=False)

  	-b reps           Specify number of bootstraps default=1000

  	-z zeta           Specify family support value, range [0, 1]

  	-m rsMethod       Specify re-sampling method: b (Bootstrap), j (Jackknife),
                    c (Permute species for each character), o (Permute
                    character order), s (Permute within species), r (Rewrite
                    data)) (default='d')

  	-s seed           Specify initial seed. if both initialSeed and loadseed
                    option are not provided then system time will be taken as
                    the default seed.

  	-l loadSeed       Specify path to a file containing previous state
                    (default=None)

  	-o outFile        Specify the file to output the results. (default=
                    Results.gfi)

  	-od outDir        Specify path to directory to store intermediate results.
                    If not specified intermediate results will be stored in a
                    temporary directory which is deleted after program
                    execution (default=None)

  	-np NP            if -np flag is provided, then no perturbed MSA will be
                    generated. Instead Gene trees will be reconstructed from
                    Bootstraps

  	-rapidnj RAPIDNJ  if -rapidnj flag is provided, then Gene trees will be
                    reconstructed from Bootstraps using RapidNJ tool

  	-fastest FASTEST  speed up the neighbor joining phase & reduce memory usage
                    (recommended for >50,000 sequences)

  	-gtr GTR          generalized time-reversible model (nucleotide alignments
                    only)

  	-wag WAG          Whelan-And-Goldman 2001 model (amino acid alignments only)

  	-gamma GAMMA      after optimizing the tree under the CAT approximation,
                    rescale the lengths to optimize the Gamma20 likelihood

