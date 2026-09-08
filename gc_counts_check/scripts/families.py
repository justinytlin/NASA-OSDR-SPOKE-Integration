"""Anatomy node families (SPOKE v2 / UBERON names) used by the checks: name -> (include regex, exclude regex)."""
FAMILIES = {
    "hippocampus": (r"hippocamp|dentate gyrus|Ammon|subiculum", None),
    "thymus": (r"thym", None),
    "spleen": (r"spleen|splenic", r"flexure"),
    "liver": (r"liver|hepat", r"hepatopancreas|flexure|arthropod|hepatic cecum"),
    # post hoc, narrower: parenchymal liver terms only (no ducts, ligaments, vessels)
    "liver parenchyma": (r"^liver$|liver lobule|hepatic cord|hepatic acinus|perisinusoidal|hepatic sinusoid|"
                         r"liver parenchyma|lobe of liver|liver (?:left|right|caudate|quadrate)|hepatocyte|central vein of liver", None),
    "cerebral cortex": (r"cerebral cortex|neocortex|isocortex|prefrontal|frontal cortex|somatosensory|visual cortex", None),
    "heart": (r"heart|cardiac|myocard", r"cardiac stomach|cardiac gastric|corpus cardiacum|lymph heart|dorsal vessel"),
    "kidney": (r"kidney|renal|nephron|glomerul", r"adrenal|suprarenal|interrenal|nephridi|Hatschek|lieno|olfactory|cerebell|preglomerular"),
    "lung": (r"\blung|pulmonary|bronch|alveolus of lung|alveolar (?:sac|wall|system|duct)", r"pulmonum"),
    "skeletal muscle": (r"skeletal muscle|gastrocnemius|quadriceps|soleus|tibialis", None),
    "bone marrow": (r"bone marrow", None),
    "skin": (r"\bskin|epidermis|dermis", r"dermatome|hypodermis"),
    "intestine": (r"intestin|duoden|jejun|ileum|colon\b|colonic|rectum", None),
    "testis": (r"\btestis|testic|seminiferous", None),
    "blood": (r"^blood$|blood plasma|peripheral blood", None),
    # added for the RR-3 payload
    "adrenal gland": (r"adrenal|suprarenal|interrenal", None),
    "eye": (r"\beye\b|eyeball|ocular|cornea|\blens\b|retina|sclera|conjunctiv|uvea|\biris\b|ciliary body|choroid|vitreous|aqueous humor",
            r"choroid plexus|choroidal artery|anterior choroidal|eyelid|eyebrow|lens placode|eye primordium|compound eye|ocellus|eyespot|"
            r"preoptic|supraoptic|optic (?:tract|radiation|tectum|lobe|chiasm)|sensillum|extra-ocular|ocular pre-muscle|eye skin|"
            r"eye sebaceous|sclerotome|scleral ossicle|retinaculum|macula densa|fovea capitis"),
    "retina": (r"retina|photoreceptor|optic disc|fovea centralis|foveola|macula lutea|parafovea|perifovea", r"retinaculum"),
    "brain": (r"brain|cerebr|cerebell|hippocamp|thalam|amygdala|striatum|telenceph|dienceph|midbrain|hindbrain|\bpons\b|medulla oblongata|"
              r"olfactory bulb|forebrain|neocortex|isocortex|pallium|basal ganglia|substantia nigra|corpus callosum",
              r"cerebral vein|cerebral artery|braincase|Myxin|ascidian|protocerebrum|cerebral ganglion"),
}
# the fixed competition set used by check_anatomy.py (the dataset's own family is added if not in it)
CORE_ORGAN_FAMILIES = ["hippocampus", "thymus", "spleen", "liver", "cerebral cortex", "heart", "kidney", "lung",
                       "skeletal muscle", "bone marrow", "skin", "intestine", "testis", "blood"]
SOURCE_FAMILY = {"hippocampus": "hippocampus", "thymus": "thymus", "spleen": "spleen", "liver": "liver"}


def family_masks(names):
    names = names.astype(str)
    out = {}
    for fam, (inc, exc) in FAMILIES.items():
        m = names.str.contains(inc, case=False, regex=True)
        if exc:
            m &= ~names.str.contains(exc, case=False, regex=True)
        out[fam] = m.values
    return out
