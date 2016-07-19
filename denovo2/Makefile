all:
	cd utils; gmake
	cd asid.1.0; gmake
	chmod 777 assemble.pl; chmod 777 analyze.pl; chmod 777 asid.1.0/asid_assembly.sh; chmod 777 asid.1.0/asid_assembly_one.sh; chmod 777 saet.2.2/saet_mp
	cd velvet_0.7.55; gmake color
	cd MUMmer3.22; gmake

clean: 
	cd utils; gmake clean
	cd asid.1.0; gmake clean
	cd velvet_0.7.55; gmake clean
	cd MUMmer3.22; gmake clean

