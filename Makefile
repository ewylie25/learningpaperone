all: b01.svg b02.svg b03.svg b04.svg b05.svg b06.svg b07.svg b08.svg b09.svg b10.svg b11.svg b12.svg b13.svg b14.svg b15.svg b16.svg b17.svg b18.svg b19.svg b20.svg b21.svg#alpha_cyclodextrin.svg cyclodextrin.svg  adeddu.svg adeddu2.svg wildsugar.svg porphyrin.svg \
	maloush.svg heme.svg heptane.svg random_db.svg random_db2.svg adeddu3.svg adeddu4.svg tropinone.svg adeddu_reaction.svg \
	adeddu_reaction1.svg adeddu_reaction2.svg adeddu5.svg #Maitotoxin.svg sucrose.svg 

smarts_stats.txt: source_SMILES/10000_smiles.txt smarts.txt
	./common_fragments.py $^ $@

random_db="Clc1cccc2cc3ccccc3cc12"
random_db.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds_to_json.py $< 10000 $(random_db) | xargs java -jar ./colorer.jar $@

random_db2="FC(F)(F)COC1C2CC3CC(C2)CC1C3"
random_db2.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(random_db2) | xargs java -jar ./colorer.jar $@

alpha_cyclodextrin="C([C@@H]1[C@@H]2[C@@H]([C@H]([C@H](O1)O[C@@H]3[C@H](O[C@@H]([C@@H]([C@H]3O)O)O[C@@H]4[C@H](O[C@@H]([C@@H]([C@H]4O)O)O[C@@H]5[C@H](O[C@@H]([C@@H]([C@H]5O)O)O[C@@H]6[C@H](O[C@@H]([C@@H]([C@H]6O)O)O[C@@H]7[C@H](O[C@H](O2)[C@@H]([C@H]7O)O)CO)CO)CO)CO)CO)O)O)O"

alpha_cyclodextrin.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(alpha_cyclodextrin) | xargs java -jar ./colorer.jar $@

maloush="O[C@@H]8[C@@H](O)[C@@H]1O[C@H](CO)[C@H]8O[C@H]7O[C@H](CO)[C@@H](O[C@H]6O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O[C@H]4O[C@H](CO)[C@@H](O[C@H]3O[C@H](CS)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O1)[C@H](O)[C@H]2O)[C@H](O)[C@H]3O)[C@H](O)[C@H]4O)[C@H](O)[C@H]5O)[C@H](O)[C@H]6O)[C@H](O)[C@H]7O"

maloush.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(maloush) | xargs java -jar ./colorer.jar $@

cyclodextrin="O[C@@H]8[C@@H](O)[C@@H]1O[C@H](CO)[C@H]8O[C@H]7O[C@H](CO)[C@@H](O[C@H]6O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O[C@H]4O[C@H](CO)[C@@H](O[C@H]3O[C@H](CO)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O1)[C@H](O)[C@H]2O)[C@H](O)[C@H]3O)[C@H](O)[C@H]4O)[C@H](O)[C@H]5O)[C@H](O)[C@H]6O)[C@H](O)[C@H]7O"

cyclodextrin.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(cyclodextrin) | xargs java -jar ./colorer.jar $@

Maitotoxin="C[C@H](CC[C@@H]([C@@H]([C@H](C)C[C@H](C(=C)/C(=C/CO)/C)O)O)OS(=O)(=O)[O-])[C@H]([C@@H](C)[C@H]1[C@@H]([C@@H]([C@H]2[C@H](O1)[C@@H](C[C@]3([C@H](O2)C[C@H]4[C@H](O3)C[C@]5([C@H](O4)[C@H]([C@H]6[C@H](O5)C[C@H]([C@H](O6)[C@@H]([C@H](C[C@H]7[C@@H]([C@@H]([C@H]8[C@H](O7)C[C@H]9[C@H](O8)C[C@H]1[C@H](O9)[C@H]([C@@H]2[C@@H](O1)[C@@H]([C@H]([C@@H](O2)[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@@H](C[C@H](C[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@H]1[C@@H](C[C@]2([C@H](O1)[C@@H]([C@]1([C@H](O2)C[C@]2([C@H](O1)CC[C@]1([C@H](O2)C[C@]2([C@H](O1)C[C@H]1[C@H](O2)CC[C@H](O1)[C@]1([C@@H](C[C@H]2[C@](O1)(C[C@H]1[C@](O2)(CC[C@]2([C@H](O1)C[C@H]1[C@](O2)(C[C@H]2[C@H](O1)C/C=C\[C@H]1[C@H](O2)C[C@H]2[C@](O1)(C[C@]1([C@H](O2)C[C@H]2[C@](O1)(CC[C@H](O2)[C@H]([C@@H](C[C@@H](C)[C@@H](C)CC=C)O)O)C)C)C)C)C)C)C)O)C)C)C)C)C)O)C)O)O)O)O)O)O)O)O)O)O)O)O)O)OS(=O)(=O)[O-])O)O)O)O)C)C)O)O)O)O.[Na+].[Na+]"

Maitotoxin.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(Maitotoxin) | xargs java -jar ./colorer.jar $@

adeddu="CCCCCCCCCCCOCOCOCOCOCOCOCOCO"
adeddu.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu) | xargs java -jar ./colorer.jar $@

adeddu2="OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCCCCCCCCCCCCCCSP(O)(O)=O"
adeddu2.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu2) | xargs java -jar ./colorer.jar $@

adeddu_reaction="Cc1ccccc1"
adeddu_reaction.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu_reaction) | xargs java -jar ./colorer.jar $@

adeddu_reaction1="Cc1cc(Cl)ccc1"
adeddu_reaction1.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu_reaction1) | xargs java -jar ./colorer.jar $@

adeddu_reaction2="Cc1cc(Cl)cc(Cl)c1"
adeddu_reaction2.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu_reaction2) | xargs java -jar ./colorer.jar $@

sucrose="O1[C@H](CO)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O[C@@]2(O[C@@H]([C@@H](O)[C@@H]2O)CO)CO"
sucrose.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(sucrose) | xargs java -jar ./colorer.jar $@

wildsugar="O=C(O[C@H]2O[C@@H]([C@@H](O[C@@H]1O[C@H](COC(=O)C)[C@@H](OC(=O)C)[C@H](OC(=O)C)[C@H]1OC(=O)C)[C@H](OC(=O)C)[C@H]2OC(=O)C)COC(=O)C)C"
wildsugar.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(wildsugar) | xargs java -jar ./colorer.jar $@

PORPHYRIN="c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"
porphyrin.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(PORPHYRIN) | xargs java -jar ./colorer.jar $@

heme="Cc1/c/2c/c3n/c(c\c4c(c(c([n-]4)/cc/5\nc(/cc(/c1CCC(=O)O)\[n-]2)C(=C5C)CCC(=O)O)C)C=C)/C(=C3C=C)C.[Fe+2]"
heme.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(heme) | xargs java -jar ./colorer.jar $@

heptane="C1CCCCCC1"
heptane.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(heptane) | xargs java -jar ./colorer.jar $@

adeddu3="SCCOCCOCCOCCOCCO"
adeddu3.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu3) | xargs java -jar ./colorer.jar $@

tropinone="CN1[C@@H]2CC[C@H]1CC(=O)C2"
tropinone.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(tropinone) | xargs java -jar ./colorer.jar $@

adeddu4="CCCCCCCCCCOCCOCCOCCOCCOCCOCCO"
adeddu4.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu4) | xargs java -jar ./colorer.jar $@

adeddu5="[CH3:1][O:2][CH2:3][CH2:4][CH2:5][CH:3]([CH2:4][CH2:5][CH:3]([CH2:4][CH2:5][CH:3]([CH2:4][CH2:5][CH:3]([CH2:4][CH2:5][CH:3]([O:2][CH3:1])[CH:4]=[CH2:5])[O:2][CH3:1])[O:2][CH3:1])[O:2][CH3:1])[O:2][CH3:1]"
adeddu5.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(adeddu5) | xargs java -jar ./colorer.jar $@
	
	
b01="CC(=CCC/C(=C/C=C/C(=C/C=C/C(=C/C=C(/C=C/C=C(/C=C/C=C(/CCC=C(C)C)\C)\C)\C)/C)/C)/C)C"
b01.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b01) | xargs java -jar ./colorer.jar $@
	
b02="C(=C\C(=C\C=C\C(=C\C=C\C=C(\C=C\C=C(\C=C\C1=C(CC(CC1(C)C)O)C)/C)/C)\C)\C)/C1=C(CC(CC1(C)C)O)C"
b02.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b02) | xargs java -jar ./colorer.jar $@

b03="P(=O)([O-])([O-])OP(=O)(OCCC(=C)C)[O-]"
b03.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b03) | xargs java -jar ./colorer.jar $@
	
b04="O=C/C=C(/CCC=C(C)C)\C"
b04.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b04) | xargs java -jar ./colorer.jar $@
	
b05="c1(c(c(c2c(c1C=O)c(c(c(c2)C)c1c(c2c(cc1C)c(c(c(c2C=O)O)O)C(C)C)O)O)C(C)C)O)O"
b05.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b05) | xargs java -jar ./colorer.jar $@
	
b06="c1(C(c2ccc(cc2)c2c(c(c(c(c2)c2ccccc2)c2ccccc2)c2ccccc2)c2ccccc2)(c2ccc(cc2)c2c(c(c(c(c2)c2ccccc2)c2ccccc2)c2ccccc2)c2ccccc2)c2ccc(cc2)c2c(c(c(c(c2)c2ccccc2)c2ccccc2)c2ccccc2)c2ccccc2)ccc(cc1)c1c(c(c(c(c1)c1ccccc1)c1ccccc1)c1ccccc1)c1ccccc1"
b06.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b06) | xargs java -jar ./colorer.jar $@
	
b07="C(=O)(CCC(=O)OC(CO)CO)OCC(OC(=O)CCC(=O)OCC(OC(=O)CCCCC(=O)OCC(OC(=O)CCCCC(=O)OC(COC(=O)CCCCC(=O)OC(COC(=O)CCC(=O)OC(COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCCCC(=O)OC(COC(=O)CCC(=O)OC(COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCCCC(=O)OC(COC(=O)CCC(=O)OC(COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO)COC(=O)CCC(=O)OC(CO)CO"
b07.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b07) | xargs java -jar ./colorer.jar $@
	
b08="c12c(cc(c(c1)C(c1c(cc(c(c1)C(c1c(cc(c(c1)C(c1cc(c(cc1O)O)C2CCCCCCCCCCC)CCCCCCCCCCC)O)O)CCCCCCCCCCC)O)O)CCCCCCCCCCC)O)O"
b08.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b08) | xargs java -jar ./colorer.jar $@
	
b09="c12c(cc(c(c1)C(c1c(cc(c(c1)C(c1c(cc(c(c1)C(c1cc(c(cc1O)O)C2OCCOCCOCCOC)OCCOCCOCCOC)OCCCCCCCCCCC)OCCCCCCCCCCC)OCCOCCOCCOC)O)O)OCCOCCOCCOC)OCCCCCCCCCCC)OCCCCCCCCCCC"
b09.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b09) | xargs java -jar ./colorer.jar $@
	
b10="c12cc(cc(c1O)Cc1cc(cc(c1O)Cc1cc(cc(c1O)Cc1c(c(cc(c1)C(C)(C)C)C2)O)C(C)(C)C)C(C)(C)C)C(C)(C)C"
b10.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b10) | xargs java -jar ./colorer.jar $@
	
b11="c12cc(cc3c1OP(=O)(Oc1c(cc(cc1C2)C(C)(C)C)Cc1cc(cc2c1OP(=O)(Oc1c(cc(cc1Cc1cc(cc4c1OP(=O)(Oc1c(cc(cc1C4)C(C)(C)C)Cc1cc(cc4c1OP(=O)(Oc1c(cc(cc1C4)C(C)(C)C)C3)Cl)C(C)(C)C)Cl)C(C)(C)C)C(C)(C)C)C2)Cl)C(C)(C)C)Cl)C(C)(C)C"
b11.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b11) | xargs java -jar ./colorer.jar $@
	
b12="C1CN2CCOCCOCCOCCN(CCOCCO1)CCCC2"
b12.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b12) | xargs java -jar ./colorer.jar $@
	
b13="C1CN2CCOCCOCCOCCN(CCOCCOCCO1)CCCC2"
b13.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b13) | xargs java -jar ./colorer.jar $@
	
b14="C1CN2CCOCCNCCOCCN(CCOCCNCCO1)CCOCCC2"
b14.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b14) | xargs java -jar ./colorer.jar $@
	
b15="C1OCCOCCOCCCOCCOCC1"
b15.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b15) | xargs java -jar ./colorer.jar $@
	
b16="[nH]1c2c(c(c1/C=C/1\N=C(/C=c/3\[nH]/c(=C\C4=N/C(=C\2)/C(=C4C=C)C)/c(c3CCC(=O)O)C)C(=C1C)CCC(=O)O)C)C=C"
b16.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b16) | xargs java -jar ./colorer.jar $@
	
b17="C(C(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
b17.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b17) | xargs java -jar ./colorer.jar $@
	
b18="c1c(c(cc(c1)/C=C/C(=O)CC(=O)/C=C/c1cc(c(cc1)O)OC)OC)O"
b18.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b18) | xargs java -jar ./colorer.jar $@
	
b19="n1c(nc(c2c1c(nc(n2)N(CCO)CCO)N1CCCCC1)N1CCCCC1)N(CCO)CCO"
b19.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b19) | xargs java -jar ./colorer.jar $@
	
b20="c1(c(cc2c(c1)Oc1c(O2)cc(c(c1)Cl)Cl)Cl)Cl"
b20.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b20) | xargs java -jar ./colorer.jar $@
	
b21="c1cc(ccc1C#Cc1c(cc(cc1C#Cc1ccc(cc1)SC)C#N)C#Cc1ccc(cc1)SC)N(c1ccc(cc1)C#Cc1c(cc(cc1C#Cc1ccc(cc1)SC)C#N)C#Cc1ccc(cc1)SC)c1ccc(cc1)C#Cc1c(cc(cc1C#Cc1ccc(cc1)SC)C#N)C#Cc1ccc(cc1)SC"
b22.svg: smarts_stats.txt weight_bonds.py
	./weight_bonds.py $< 10000 $(b22) | xargs java -jar ./colorer.jar $@
	
	

clean:
	rm -rf smarts_stats.txt *.svg
