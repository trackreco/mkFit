# Description of implementation of multiple mkFit iterations

- The branch is up-to-date with respect to devel.

- The main changes in this branch affect the following files:

(1) mkFit/SteeringParams.h:
 - three additional classes, which are iteration-dependent:
 (a) IterationParams: a container for 'nlayers_per_seed', 'maxCandsPerSeed', 'maxHolesPerCand', 'maxConsecHoles', and 'chi2Cut';
 (b) IterationLayerConfig: a container for layer-specific iteration-dependent configurations (e.g., hit selection windows);
 (c) IterationConfig: a container for all of the above, including a virtual functions to import seeds (import_seeds)
 - one additional struct, which is iteration-dependent:
 (a) MkSeedPacket: a container of iteration-specific event objects ('m_seedEtaSeparators_', 'm_seedMinLastLayer_', 'm_seedMaxLastLayer_', 'm_layerHits_', 'm_inseeds', 'm_outtrks')

(2) Geoms/CMS-2017.cc:
 - an instance of IterationConfig is created in Geoms/CMS-2017.cc, to be passed to MkBuilder constructor, which sets all iteration-dependent objects/parameters.

(3) mkFit/MkBuilder[.cc,.h]
 - all iteration-dependent parameters (regions, steering parameters, etc.) are moved out of MkBuilder, and into IterationConfig, which must be passed to MkBuilder constructor to have one MkBuilder per iteration.
