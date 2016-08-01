###################
## Cluster stuff ##
###################

SUBMIT = eval

#################
## Directories ##
#################

SRCDIR = $(CURDIR)/src
INPUT = $(CURDIR)/input
MUTATION = $(CURDIR)/mutations
PLOTDATA = $(CURDIR)/plotdata
PLOTDIR = $(CURDIR)/plots
FIGUREDIR = $(CURDIR)/figures

############
## Files ##
############

GENOME = genome.gbk

UNIPROT = $(MUTATION)/uniprot.txt
UNIPROTSIZES = $(MUTATION)/uniprot_sizes.txt

# Uniprot features
FEATURES = $(MUTATION)/features.tsv
FEATURECATEGS = $(INPUT)/features.txt
FILTEREDFEATURES = $(MUTATION)/filtered_features.tsv
FEATURESBED = $(MUTATION)/features.bed

# Uniprot mutations
MUTATIONS = $(MUTATION)/mutations.txt
ANNMUTATIONS = $(INPUT)/mutations.tsv
TOLMUTATIONS = $(MUTATION)/tolerated.txt
DELMUTATIONS = $(MUTATION)/deleterious.txt

# All nonsyn mutations observed (external data)
ALLNONSYN = $(INPUT)/all_nonsynmuts.tsv
ALLNONSYNBED = $(MUTATION)/all_nonsynmuts.bed
# All sift scores for the mutations
ALLSIFT = $(INPUT)/all_sift.tsv
ALLSIFTBED = $(MUTATION)/all_sift.bed
OBSSIFT = $(INPUT)/all_sift_nonsyn.tsv
OBSSIFTBED = $(MUTATION)/all_sift_nonsyn.bed
# All foldx scores for the mutations
ALLFOLDX = $(INPUT)/all_foldx.tsv
ALLFOLDXBED = $(MUTATION)/all_foldx.bed
OBSFOLDX = $(INPUT)/all_foldx_nonsyn.tsv
OBSFOLDXBED = $(MUTATION)/all_foldx_nonsyn.bed
# Accessibility data
ACCESSIBILITY = $(INPUT)/all_accessibility.bed
ALLACCESSIBILITY = $(MUTATION)/all_accessibility.bed

# Overlap between Uniprot features and non-syn mutations
OBSFEATURES = $(MUTATION)/observed_features.bed
OBSOTHERS = $(MUTATION)/observed_others.bed
# SIFT data that overlaps with Uniprot features
ALLSIFTFEATURES = $(MUTATION)/all_sift_features.bed
ALLSIFTOTHERS = $(MUTATION)/all_sift_others.bed
SIFTFEATURES = $(MUTATION)/sift_features.bed
SIFTOTHERS = $(MUTATION)/sift_others.bed

ESSENTIAL = $(MUTATION)/essential.txt
ALLESSENTIALMUTS = $(MUTATION)/essential_nonsynmuts.tsv

ECKFILE = $(MUTATION)/eck_uniprot.tsv
GIFILE = $(MUTATION)/gi_uniprot.tsv

CONVERSION = $(MUTATION)/locus_uniprot.tsv

TOLMUTATIONS = $(MUTATION)/tolerated.txt
DELMUTATIONS = $(MUTATION)/deleterious.txt
# TODO: add rule to generate this from observed SIFT data
TOLSIFT = $(INPUT)/tolerated.sift.txt
DELSIFT = $(INPUT)/deleterious.sift.txt

TREE = $(INPUT)/tree.nwk
NONSYNCOUNT = $(INPUT)/strains_nonsyn.txt
PANGENOMECOUNT = $(INPUT)/strains_pangenome.txt
EVOLUTION = $(INPUT)/evolution_experiment.txt

FEATURESDATA = $(PLOTDATA)/features.tsv
SIFTFEATURESDATA = $(PLOTDATA)/sift_features.tsv
ACCESSIBILITYDATA = $(PLOTDATA)/accessibility.tsv
FOLDXACCESSIBILITYDATA = $(PLOTDATA)/foldx_accessibility.tsv

TREEPLOT = $(PLOTDIR)/tree.svg
TREEBARS = $(PLOTDIR)/tree_colorbars.svg
TREELEGEND = $(PLOTDIR)/tree_legend.svg
CONSTRAINTSPLOT = $(PLOTDIR)/constraints.svg

FIGURE1 = $(FIGUREDIR)/figure_1.svg

##############################
## Non-synonymous mutations ##
##############################

$(ALLNONSYNBED): $(ALLNONSYN)
	$(SRCDIR)/nonsyn2bed $(ALLNONSYN) > $@

##################
## Uniprot data ##
##################

$(UNIPROT): $(GENOME)
	$(SRCDIR)/gbk2uniprot $(GENOME) > $@
$(UNIPROTSIZES): $(GENOME)
	$(SRCDIR)/gbk2uniprot $(GENOME) --size > $@

$(FEATURES): $(UNIPROT)
	$(SRCDIR)/uniprot_features $< > $@

$(FILTEREDFEATURES): $(FEATURES) $(FEATURECATEGS)
	-rm $@
	for categ in $$(cat $(FEATURECATEGS)); do grep $$categ $(FEATURES) >> $@; done

$(FEATURESBED): $(FILTEREDFEATURES)
	$(SRCDIR)/features2bed $(FILTEREDFEATURES) > $(FEATURESBED)

$(OBSFEATURES): $(ALLNONSYNBED) $(FEATURESBED)
	bedtools intersect -a $(ALLNONSYNBED) -b $(FEATURESBED) | sort > $@
$(OBSOTHERS): $(ALLNONSYNBED) $(FEATURESBED)
	bedtools intersect -a $(ALLNONSYNBED) -b $(FEATURESBED) -v | sort > $@

$(MUTATIONS): $(UNIPROT)
	$(SRCDIR)/uniprot_mutations $(UNIPROT) > $(MUTATIONS)

$(TOLMUTATIONS): $(ANNMUTATIONS) $(ALLESSENTIALMUTS)
	$(SRCDIR)/mutations2nonsyn $(ANNMUTATIONS) $(TOLMUTATIONS) $(DELMUTATIONS) --essential $(ALLESSENTIALMUTS) 

###############
## SIFT data ##
###############

$(ALLSIFTBED): $(ALLSIFT)
	$(SRCDIR)/sift2bed $< > $@
$(OBSSIFTBED): $(OBSSIFT)
	$(SRCDIR)/sift2bed $< > $@

$(ALLSIFTFEATURES): $(FEATURESBED) $(ALLSIFTBED)
	bedtools intersect -a $(ALLSIFTBED) -b $(FEATURESBED) > $@
$(ALLSIFTOTHERS): $(FEATURESBED) $(ALLSIFTBED)
	bedtools intersect -a $(ALLSIFTBED) -b $(FEATURESBED) -v > $@
$(SIFTFEATURES): $(FEATURESBED) $(OBSSIFTBED)
	bedtools intersect -a $(OBSSIFTBED) -b $(FEATURESBED) > $@
$(SIFTOTHERS): $(FEATURESBED) $(OBSSIFTBED)
	bedtools intersect -a $(OBSSIFTBED) -b $(FEATURESBED) -v > $@

################
## FOLDX data ##
################

$(ALLFOLDXBED): $(ALLFOLDX)
	$(SRCDIR)/sift2bed $< > $@
$(OBSFOLDXBED): $(OBSFOLDX)
	$(SRCDIR)/sift2bed $< > $@

########################
## Accessibility data ##
########################

$(ALLACCESSIBILITY): $(ACCESSIBILITY) $(ALLNONSYNBED) 
	bedtools intersect -a $(ACCESSIBILITY) -b $(ALLNONSYNBED) > $@

##########################
## Locus to Uniprot IDs ##
##########################

$(CONVERSION): $(GENOME)
	$(SRCDIR)/gbk2locusuniprot $(GENOME) > $(CONVERSION)

########################
## ECK TO UNIPROT IDS ##
########################

$(ECKFILE): $(GENOME)
	$(SRCDIR)/eck2uniprot $(GENOME) | sort | uniq > $(ECKFILE)

#######################
## GI TO UNIPROT IDS ##
#######################

$(GIFILE): $(GENOME)
	$(SRCDIR)/gbk2giuniprot $(GENOME) | sort | uniq > $(GIFILE)

#######################
## Gene essentiality ##
#######################

$(ESSENTIAL): $(CONVERSION)
	wget -O $(MUTATION)/PECData.dat "http://shigen.nig.ac.jp/ecoli/pec/download/files/PECData.dat" && \
	awk -F'\t' '{if ($$10 == 1) print $$4"\t"$$10}' $(MUTATION)/PECData.dat | awk -F',' '{print $$1}' | sort | uniq > $(MUTATION)/PECData.txt && \
	$(SRCDIR)/essential2uniprot $(MUTATION)/PECData.txt $(CONVERSION) > $@

$(ALLESSENTIALMUTS): $(ALLNONSYN) $(ESSENTIAL)
	-rm $(ALLESSENTIALMUTS)
	-for essential in $$(cat $(ESSENTIAL)); do grep $$essential $(ALLNONSYN) >> $(ALLESSENTIALMUTS); done

##########################
## Plot data generation ##
##########################

$(FEATURESDATA): $(ESSENTIAL) $(UNIPROTSIZES) $(FEATURESBED) $(OBSFEATURES) $(OBSOTHERS)
	$(SRCDIR)/run_constraints_features $(ESSENTIAL) $(UNIPROTSIZES) $(FEATURESBED) $(OBSFEATURES) $(OBSOTHERS) --bootstraps 100 > $@

$(SIFTFEATURESDATA): $(ESSENTIAL) $(SIFTFEATURES) $(SIFTOTHERS) $(ALLSIFTFEATURES) $(ALLSIFTOTHERS)
	$(SRCDIR)/run_constraints_features_sift $(ESSENTIAL)  $(SIFTFEATURES) $(SIFTOTHERS) $(ALLSIFTFEATURES) $(ALLSIFTOTHERS) --bootstraps 100 > $@

$(ACCESSIBILITYDATA): $(ESSENTIAL) $(ACCESSIBILITY) $(ALLACCESSIBILITY)
	$(SRCDIR)/run_constraints_accessibility $(ESSENTIAL) $(ACCESSIBILITY) $(ALLACCESSIBILITY) --bootstraps 100 --bins 100 > $@

$(FOLDXACCESSIBILITYDATA): $(ESSENTIAL) $(ACCESSIBILITY) $(ALLFOLDXBED) $(OBSFOLDXBED)
	$(SRCDIR)/run_constraints_accessibility_foldx $(ESSENTIAL) $(ACCESSIBILITY) $(ALLFOLDXBED) $(OBSFOLDXBED) --bootstraps 100 --buried 50 > $@

######################
## Plots generation ##
######################

$(TREEPLOT): $(TREE) $(EVOLUTION) $(NONSYNCOUNT) $(PANGENOMECOUNT)
	$(SRCDIR)/run_tree_generation $(TREE) $(EVOLUTION) $(NONSYNCOUNT) $(PANGENOMECOUNT) $@ --height 7 --dpi 90

$(TREEBARS): $(TREE) $(NONSYNCOUNT) $(PANGENOMECOUNT)
	$(SRCDIR)/run_tree_colorbar $(TREE) $(NONSYNCOUNT) $(PANGENOMECOUNT) $@ --height 0.3 --width 5.7 --dpi 90

$(TREELEGEND):
	$(SRCDIR)/run_tree_legend $@ --height 0.5 --width 2 --dpi 90

$(CONSTRAINTSPLOT): $(FEATURESDATA) $(SIFTFEATURESDATA) $(ACCESSIBILITYDATA) $(FOLDXACCESSIBILITYDATA)
	$(SRCDIR)/run_constraints_plot $(FEATURESDATA) $(SIFTFEATURESDATA) $(ACCESSIBILITYDATA) $(FOLDXACCESSIBILITYDATA) $@ --height 7 --width 7 --dpi 90

###############
## Figures 1 ##
###############

$(FIGURE1): $(TREEPLOT) $(TREEBARS) $(TREELEGEND) $(CONSTRAINTSPLOT)
	$(SRCDIR)/run_figure_1 $(TREEPLOT) $(TREEBARS) $(TREELEGEND) $(CONSTRAINTSPLOT) $@

########################
## Targets definition ##
########################

constraints: $(FEATURESDATA) $(SIFTFEATURESDATA) $(ACCESSIBILITYDATA) $(FOLDXACCESSIBILITYDATA)
plots: $(TREEPLOT) $(TREEBARS) $(TREELEGEND) $(CONSTRAINTSPLOT)
figures: $(FIGURE1)

all: constraints plots figures

.PHONY: all constraints plots figures
