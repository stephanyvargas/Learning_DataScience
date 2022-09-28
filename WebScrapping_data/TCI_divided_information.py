
#Naming conventions in different platforms for the same compound
compound_identifications = ['name',  'CAS', 'code', 'grade', 'product number', 'cas rn',
              'reaxys registry number', 'pubchem substance id', 'SMILES_by_PubChem',
              'merck index (14)', 'mdl number', 'sdbs (aist spectral db)',
              'related cas rn', 'colour index', 'enzyme commission number']


#Available stock in different areas (Hyogo, Saitama and others) with price and amount
available_stock = [ 'code', '0.1ML','0.2ML','0.5ML','100G','100MG','100ML',
                    '100TUB', '10G', '10MG', '10ML', '12CAP', '12ML', '15ML',
                    '1CLM','1EACH','1G','1KG','1KIT','1MG','1ML','1SAMPLE',
                    '1SET','1VIAL','200G','200MG','200ML','20MG','20ML','250G',
                    '250MG','250ML','25G','25MG','25ML','2G','2MG','2ML',
                    '300G','300ML','3CLM','400G','400ML','500G','500MG',
                    '500ML','50G','50MG','50ML','5CLM','5EACH','5G','5MG',
                    '5ML','8ML', '90G', '96CAP' ]


#General information of the product
general_information=['code', 'Molecular Formula', 'Molecular Weight',
 'purity / analysis method', 'appearance', 'solubility in water',
 'molecular weight', 'sensitiveness', 'average active oxygen',
 'the average value of n', 'the average value of m+n',
 'the average value of (n+m)', 'content(na,drying substance)',
 'etherification value( as drying substance)', 'cw (scheduled material type)',
 'assay of mono ester', 'assay of diester']
suitability_for_tests_and_analysis = [  'suitability for absorptiometry',
 'suitability for aldehyde-analysis', 'suitability for amino acid analysis',
 'suitability for arsenic-analysis', 'suitability for berylium-analysis',
 'suitability for calcium-analysis', 'suitability for chrome-analysis',
 'suitability for cobalt-analysis', 'suitability for cyan-analysis',
 'suitability for electrophoresis', 'suitability for formaldehyde-analysys',
 'suitability for gc-analysis', 'suitability for iron-analysis',
 'suitability for lc-ms analysis', 'suitability for magnesium-analysis',
 'suitability for mass-analysis calibratio', 'suitability for melamine-analysis',
 'suitability for nitrate-analysis', 'suitability for nmr-analysis',
 'suitability for p-cresol-analysis', 'suitability for protein analysis',
 'suitability for protein-analysis', 'suitability for protein-analysis (e.coli)',
 'suitability for protein-analysis (yeast)', 'suitability for redox reagent',
 'suitability for silver-analysis', 'suitability for sulfate-analysis',
 'suitability for sulfide-analysis', 'suitability for sulfite-analysis',
 'suitability for titanium-analysis', 'suitability for use in elisa tests',
 'suitability for vanadium-analysis', 'suitability for vcm-analysis',
 'suitability test', 'suitability test for protein-analysis',
 'elemental analysis(nitrogen)', 'elemental analysis(carbon)',
 'elemental analysis(oxygen)', 'elemental analysis(sulfuride)',
 'elemental analysis(hydrogen)']
storage = [ 'store under inert gas', 'storage temperature']


#Specific details of a product
specifications_general = ['code', 'loading', 'melting range', 'drying loss', 'acid value',
 'mesomorphic range', 'viscosity (1% in water)', 'ph (1% in water)',
 'saponification value', 'degree of deuteration',
 'evaporation residue', 'specific gravity', 'sulfate',
 'homo level', 'entantiometric excess (ee)', 'molar ratio maleimido',
 'deacetylation value', 'average alkylchain length', 'limiting viscosity',
 'transmittance', 'loading (as pd)', 'electron mobility(mufet)',
 'mass concentration(au)', 'particle concentration', 'optical density']
purity = [ 'purity(nonaqueous titration)',  'purity(argentometric titration)',
 'purity(hplc)', 'purity(sodium hypochlorite method)', 'purity(gc)',
 'purity(neutralization titration)', 'purity(with total nitrogen)',
 'purity( potassium iodate method)', 'purity(periodic acid method)',
 'purity(hplc)(cad)', 'purity(bromination method)', 'purity(tpbna method)',
 'purity(chelometric titration)','purity(volhard method)', 'purity(iodometric titration)',
 'purity(volumetric analysis)', 'purity(ion exchange titration)',
 'purity(uv-vis method)', 'purity(precipitation titration)',
 'purity(methanolysis method)', 'purity(morpholine method)', 'purity(with ignition residue)',
 'purity(potassium permanganate method)', 'purity(qnmr)', 'purity(nh4scn method)',
 'purity(hplc)(ri)', 'purity(anilin method)', 'purity(iodometric back titration)',
 'purity', 'purity(formol titration)', 'purity(gasburet method)', 'purity(chelometric back titration)',
 'purity(potassium bromate method)', 'purity(ester value)', 'purity(gravimetric)',
 'purity(with total sulfur)', 'purity(gravimetric method)', 'purity(redox method)',
 'purity(neutralization back titration)', 'purity(cerium redox method)',
 'purity(butylamine method)', 'purity(nmr)', 'purity(oxime formation)',
 'purity(phase splitting method)', 'purity(saponification value)',
 'purity(trace metal basis)', 'purity(tpo method)', 'purity(1h-qnmr)',
 'purity(coupling titration)', 'purity(redox titration)', 'purity(anilide method)']
concentration = [ 'concentration', 'concentration (calcd. by abs)',
 'concentration (redox titration)', 'concentration of nah',
 'concentration( potassium iodate method)', 'concentration(argentometric back titration)',
 'concentration(argentometric titration)', 'concentration(bromination method)',
 'concentration(by measurement of drying weight)', 'concentration(chelometric method)',
 'concentration(chelometric method)unit2', 'concentration(chelometric titration)',
 'concentration(gasburet method)', 'concentration(gc)',
 'concentration(gravimetric method)', 'concentration(hnmr)',
 'concentration(iodometric back titration)', 'concentration(iodometric titration)',
 'concentration(lowry method)', 'concentration(morpholine method)',
 'concentration(nbs method)', 'concentration(neutralization back titration)',
 'concentration(neutralization titration)', 'concentration(neutralization titration) unit1',
 'concentration(nonaqueous titration)', 'concentration(oxime formation)',
 'concentration(phase splitting method)', 'concentration(potassium permanganate method)',
 'concentration(precipitation titration)', 'concentration(redox method)',
 'concentration(redox titration)', 'concentration(sec butanol method)',
 'concentration(titration)', 'concentration(with evaporation residue)',
 'concentration(with total nitrogen)']
test = [ 'functionality test', 'color test', 'effect test of tms derivatization',
 'enzyme inhibition test', 'blotting test', 'protein staining test',
 'dna staining test', 'fucose inhibition test',
 'protein stabilization test', 'performance test', 'protein labeling test(qualitative method)',
 'beta-galactosidase detection test', 'beta-glucuronidase detection test',
 'protein impurity test(sds-page)', 'chemiluminescence test (luminol/pod/h2o2)',
 'chemiluminescence test (superoxide radical)', 'conjugate test',
 'protein determination test', 'enzyme detection test',
 'ethanol insoluble matter', 'effect test of hplc derivatization',
 'high boiling impurity of silylation agent', 'effect test of gc derivalization',
 'alkaline phosphatase detection test', 'substances darkend by h2so4',
 'functionality test(tlc)', 'delayed emission','mass spectroscopy',
 'amino acid sequencing', 'interference pigment', 'inhibitor mixed test',
 'germfree test', 'titer test', 'effect test of ester derivatization',
 'phosphatase inhibition test', 'nitrite ion detection test（fluorescence method）',
 'collagen detection test', 'ir', 'alkali', 'fluorescence test',
 'aldehyde', 'free acid', 'protein denaturation test']
content_of_other_products = [ 'ammonium chloride', 'water',
 'total nitrogen', 'histidine', 'hydrogendioxide detection test(peroxidase method)',
 'sodium chloride', 'ash content', 'content of magnesium', 'sulfuric acid',
 'iron', 'heavy metals(as pb)', 'content of calcium', 'residual solvent',
 'peroxidase detection test', 'calcium', 'content (palladium)',
 'cobalt', 'ethanol', 'ammonium', 'arsenic(as)', 'total sulfur', 'content(cu)',
 'zinc', 'benzene sulfonic acid', 'acid(as pyromellitic acid)', 'triuret',
'chloride', 'boric acid', 'sodium sulfate', 'base(as koh)', 'water(vaporization method)',
 'diphenyl sulfone', 'bleomycina2', 'bleomycinb2', 'bismuth', 'toluene',
 'fluorescein', 'bis chloromethyl ether(gc)', 'heavy metals', 'phosphate',
 'dicyandiamide', 'creatinine', 'm-cresol', 'p-cresol', 'iodine value',
 'non volatile matter', 'free chlorine', 'thiamazole', 'residual chloroform',
 '3-chloro-1-butene', 'content(copper)', 'cell staining test', 'rhodium',
 'content of n-ag', 'benzene', 'iridium', 'content of cobalt', 'content(sodium)',
 'content of carbon', 'ignition residue(as sno2)', '2-methyl-2-butene',
 'free amine(as dimethylamine)', 'ammonia', 'methylamine', 'trimethylamine',
 'phenol', 'high boiling impurities', 'bromide', 'arsenic', 'content of sulfur',
 'polychlorinated biphenyl', 'ethanol(nmr)', 'chloride content', 'calcium carbonate',
 'sulfur', 'content of hf(neutralization titration)', 'sodium', 'potassium',
 'beta-1,3-1,6-glucan', 'oxalic acid', 'graphite residue', 'manganese',
 'other amino acids', 'l-glutamic acid', 'glycine', 'total acid',
 'methylene chloride(nmr)', 'hexane', 'concentration(neutralization titration) unit1',
 'haloid', 'water(value)', 'bromide content', 'iodine content', 'methanol']
fatty_acid =  ['fatty acid composition', 'stearic acid(composition of fatty acid)']
ignition_residue = [ 'ignition residue', 'ignition residue(sulfate)']
nmr = [ 'nmr(1h)', 'nmr(13c)', 'nmr']
abs_ratio = [ 'abs-ratio(a/b)', 'abs-ratio(2a/2b)', 'abs-ratio(3a/3b)',
 'abs-ratio(4a/4b)']
titer = [ 'titer', 'titer(elisa)', 'titer(fcm)']
optical_purity = [ 'optical purity', 'optical purity(lc)', 'optical purity(gc)',
 'optical purity(gc,mosher)', 'optical purity(hplc)']
Physical_state = [ 'diameter', 'length', 'median size', 'shape', 'physical state (20 deg.c)',
 'particle size(d50)', 'thick', 'average size', 'specific surface']
PH = [ 'ph', 'neutralization value']
Molar_extinction = [ 'molar extinction coefficient2', 'molar extinction coefficient',
 'molar extinction coefficient3', 'molar extinction coefficient(co complex)']


 #Specific properties of a product
properties_general = ['code', 'flash point',  'specific gravity (20/20)',
 'density(20deg.c)', 'transition interval(ph)', 'ester value',
 'boiling point', 'viscosity', 'average molecular weight',
 'freezing point', 'binding capacity', 'biotinylation ratio', 'lumo level',
 'hole mobility(mu fet)', 'exchange capacity']
absorbance_wavelength = [ 'maximum absorption wavelength', 'absorbance(275nm)',
 'absorbance(260nm)', 'absorbance(270nm)', 'absorbance(280nm)', 'absorbance(400nm)',
 'absorbance(330nm)', 'absorbance(340nm)', 'absorbance(350nm)', 'absorbance(360nm)',
 'absorbance(e1%1cm)', 'absorbance(254nm)', 'absorbance(300nm)', 'absorbance(310nm)',
 'absorbance(320nm)', 'absorbance(290nm)', 'absorbance(210nm)', 'absorbance(220nm)',
 'absorbance(230nm)', 'absorbance(e1%1cm)2', 'absorbance(370nm)', 'absorbance(380nm)',
 'absorbance(390nm)', 'absorbance(450nm)', 'absorbance(240nm)', 'absorbance(250nm)',
 'absorbance(e10%1cm)']
refractive_index = [ 'refractive index', 'refractive index n20/d']
melting_point = [ 'melting point', 'melting point(decomposition)']
solubility = [ 'degree of solubility in water', 'solubility (miscible with)',
 'solubility (insoluble in)', 'solubility (very soluble in)',
 'solubility (soluble in)', 'solubility in methanol', 'solubility (very slightly)',
 'solubility in hot etoh', 'solubility in etoh', 'solubility (slightly sol. in)',
 'solubility in hot water', 'solubility in hot 1mol/l hcl',
 'solubility in sodium hydroxide solution', 'solubility in hcl(1+1)',
 'solubility in toluene', 'solubility in hot toluene', 'solubility in hot methanol',
 'solubility in dilute hcl', 'solubility in hcl', 'solubility in hot dilute hcl',
 'solubility in hcl(1+3)', 'solubility in hot hcl(1+3)', 'solubility in 1mol/l hcl',
 'solubility in 1mol/l naoh', 'solubility in hcl(1+10)', 'solubility in acetic acid',
 'solubility in acetone', 'solubility in pyridine', 'solubility in n,n-dmf',
 'solubility in 5mol/l hcl', 'solubility in acetonitrile', 'solubility in chloroform',
 'solubility in thf', 'solubility in hot acetonitrile', 'solubility in 0.1mol/l naoh',
 'solubility in hot acetic acid', 'solubility in hot pyridine', 'solubility in hot etoh(50vol%)',
 'solubility in na2co3', 'solubility in ethylacetate', 'solubility in hot acetone',
 'solubility in hot dioxane', 'solubility in hot 1mol/l naoh', 'solubility in hot hcl',
 'solubility in hot 0.1m hcl', 'solubility in etoh(95vol%)', 'solubility in 0.5mol/l hydrochloric acid',
 'solubility in dioxane', 'solubility in etoh(50vol%)', 'solubility in hot dmf',
 'solubility in hot mek', 'solubility in dichloromethane', 'solubility in nh3ap.(2+3)',
 'solubility in hot chloroform', 'solubility in h2so4(1+1)', 'solubility in hot thf',
 'solubility in nh3aq.', 'solubility in formic acid', 'solubility in naoh(100g/l)',
 'solubility in 0.2mol/l naoh', 'solubility in 1-methyl-2-pyrrolidone',
 'solubility in toluene-etoh mix.', 'solubility in 2-propanol']
specific_rotation =   ['specific rotation', 'specific rotation [a]20/d',
 'specific rotation(value)', 'specific rotation [a]25/d', 'specific rotation [a]25']
absorbance = [ 'absorbance', 'absorbance of cu complex',
 'molar absorbance(al complex)', 'absorbance2']
lambda_max = [ 'lambda max.', 'lambda max.2', 'lambda max.1', 'lambda max.3']


#Transportation information
Transport_Information = ['code', 'Shipment Information', 'un number',
 'Packaging and container', 'packing group', 'class', 'air transportation']


#Categories were the product is more likely to be used
Product_categories = ['code', 'Life Science', 'Chemistry', 'Chemicals by Class',
 'Materials Science', 'Analytical Chemistry']


 #Safety and precaution guidelines
GHS_precautionary_statement = ['code', 'condition to avoid',
 'signal word', 'poisonous or deleterious']


 #Laws that are governed by agreeing to handle this product
Related_Laws = ['code', 'Chemical Substance Law_No', 'rtecs#',
 'fire defense law', 'prtr law (new specific chemical)',
 'narcotics and psychotropics control law', 'ishl',
 'chemical substance law (encs)', 'pharmaceutical affairs law (scheduled)',
 'protection of the ozone layer law (type (specified material))']


 #Entries that quite did not fit other tables
Other = ['active alkali(libu)', 'free p-nitrophenol',
 'nitorogen of cyanamide', 'oxirane oxygen', 'optical isomer',
 'petroleum ether soluble matter', 'disodium alpha-glycerophosphate',
 'reducing suger', 'human serum albmin (hsa) binding activity',
 'competitive avtivity in competitive drug',
 'non-competitive avtivity in non-competitive drug', 'se']
