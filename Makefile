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

# Overlap between Uniprot features and non-syn mutations
OBSFEATURES = $(MUTATION)/observed_features.bed
OBSOTHERS = $(MUTATION)/observed_others.bed
# SIFT data that overlaps with Uniprot features
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
TOLSIFT = $(MUTATION)/tolerated.sift.txt
DELSIFT = $(MUTATION)/deleterious.sift.txt

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

$(SIFTFEATURES): $(FEATURESBED) $(ALLSIFTBED)
	bedtools intersect -a $(ALLSIFTBED) -b $(FEATURESBED) > $@
$(SIFTOTHERS): $(FEATURESBED) $(ALLSIFTBED)
	bedtools intersect -a $(ALLSIFTBED) -b $(FEATURESBED) -v > $@

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
	curl http://ogeedb.embl.de/downloads/83333_dataset353.txt.gz > $(MUTATION)/83333_dataset353.txt.gz && \
	gunzip -f $(MUTATION)/83333_dataset353.txt.gz && \
	$(SRCDIR)/essential2uniprot $(MUTATION)/83333_dataset353.txt $(CONVERSION) > $@

$(ALLESSENTIALMUTS): $(ALLNONSYN) $(ESSENTIAL)
	-rm $(ALLESSENTIALMUTS)
	for essential in $$(cat $(ESSENTIAL)); do grep $$essential $(ALLNONSYN) >> $(ALLESSENTIALMUTS); done

########################
## Targets definition ##
########################

features: $(FEATURES)
mutations: $(MUTATIONS)
annmutations: $(TOLMUTATIONS)
test: $(UNIPROTSIZES) $(OBSFEATURES) $(OBSOTHERS) $(SIFTFEATURES) $(SIFTOTHERS)

all: features mutations annmutations

.PHONY: all features mutations annmutations test
