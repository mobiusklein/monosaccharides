import io
import re
import datetime
import zlib
import hashlib

import glypy
from glypy.structure.glycan_composition import from_iupac_lite, to_iupac_lite
from glypy.structure.stereochemistry import stereocode
from glypy.io import glycoct, iupac, wurcs
from glypy.composition import formula
from glypy.io.nomenclature import identity, synonyms

from psims.controlled_vocabulary import entity, obo, ControlledVocabulary, obj_to_xsdtype

parser = obo.OBOParser(io.BytesIO())

seen = set()
for lab, scls in glypy.structure.SuperClass:
    if scls.value < 11 and scls not in seen:
        seen.add(scls)

supercls_entities = sorted(seen)
supercls_to_id = {}
for t in supercls_entities:
    if t.name != 'x':
        mono = (from_iupac_lite(t.name.title()))
        mono_tp = {
            "id": "MONO:%s" % hex(abs(zlib.crc32(str(mono) + formula(mono.total_composition())))).upper()[2:].zfill(8),
            "name": str(mono),
            "def": "A generic monosaccharide with %d backbone carbons" % (mono.superclass.value),
            "synonyms": [
                dialect(mono) for dialect in [glycoct.dumps, iupac.dumps, wurcs.dumps,
                                              lambda x: to_iupac_lite(x)]
            ],
            "property_value": [
                "has_chemical_formula \"%s\" %s" % (formula(mono.total_composition()),
                                             obj_to_xsdtype(unicode(formula(mono.total_composition())))),
                "has_monoisotopic_mass \"%s\" %s" % (mono.mass(), obj_to_xsdtype(mono.mass())),
            ]
        }
        supercls_to_id[mono.superclass] = mono_tp['id']
        parser.current_term = mono_tp
        parser.pack()

seen = set()
handled_monos = set()
# for lab, stem in glypy.structure.Stem:
#     if stem in seen:
#         continue
#     seen.add(stem)
#     if stem.name == 'x':
#         continue
#     try:
#         mono = from_iupac_lite(stem.name.title())
#         handled_monos.add(mono)
#         mono_tp = {
#             "id": "MONO:%s" % hex(abs(zlib.crc32(str(mono) + formula(mono.total_composition())))).upper()[2:].zfill(8),
#             "name": str(mono),
#             "def": "A monosaccharide with %d backbone carbons and characteristic stereoisometry %r" % (
#                 mono.superclass.value, '-'.join(map(lambda x: str(x).title(), mono.stem))),
#             "synonyms": [
#                 dialect(mono) for dialect in [glycoct.dumps, iupac.dumps, wurcs.dumps,
#                                               lambda x: to_iupac_lite(x)]
#                 if dialect(mono) != str(mono)
#             ] + [s for s in synonyms.monosaccharides[str(mono)] if s != str(mono)],
#             "property_value": [
#                 "has_chemical_formula \"%s\" %s" % (formula(mono.total_composition()),
#                                                     obj_to_xsdtype(unicode(formula(mono.total_composition())))),
#                 "has_monoisotopic_mass \"%s\" %s" % (mono.mass(), obj_to_xsdtype(mono.mass())),
#             ]
#         }
#         parser.current_term = (mono_tp)
#         parser.pack()
#     except iupac.IUPACError as err:
#         print(err)

for label in ["dHex", "Fuc", "HexN", "HexNAc", "HexS", "HexP", "HexNAc(S)",
              "NeuAc", "NeuGc", "Neu", "HexNS", "6-aHex", "en,6-aHex"]:
    try:
        mono = from_iupac_lite(label)
#         handled_monos.add(mono)
        mono_tp = {
            "id": "MONO:%s" % hex(abs(zlib.crc32(str(mono) + formula(mono.total_composition())))).upper()[2:].zfill(8),
            "name": str(mono),
            "def": str(mono),
            "synonyms": [
                dialect(mono) for dialect in [glycoct.dumps, iupac.dumps, wurcs.dumps,
                                              lambda x: to_iupac_lite(x)]
                if dialect(mono) != str(mono)
            ] + [s for s in synonyms.monosaccharides.get(str(mono), []) if s != str(mono)] + [label],
            "property_value": [
                "has_chemical_formula \"%s\" %s" % (formula(mono.total_composition()),
                                                    obj_to_xsdtype(unicode(formula(mono.total_composition())))),
                "has_monoisotopic_mass \"%s\" %s" % (mono.mass(), obj_to_xsdtype(mono.mass())),
            ]
        }
        parser.current_term = (mono_tp)
        parser.pack()
    except iupac.IUPACError as err:
        print(err)

parser._connect_parents()
parser._simplify_header_information()

cv = ControlledVocabulary(parser.terms)

def write_header(self, header, stream):
    for key, value in header:
        stream.write("%s: %s\n" % (key, value))
    stream.write("\n")
    stream.write("\n")

def write_term(self, term, stream):
    stream.write("[Term]\nid: %s\nname: %s\ndef: \"%s\"\n" % (term.id, term.name, term.definition))
#     for xref in term.get('xref', []):
#         stream.write("xref: ")
    seen = set()
    for syn in term.get('synonyms', []):
        if syn in seen:
            continue
        seen.add(syn)
        stream.write("synonym: \"%s\" EXACT\n" % str(syn).replace("\n", "\\n"))
    for prop in term.get('property_value', []):
        stream.write("property_value: %s\n" % prop)
    stream.write("\n")

buff = io.BytesIO()

header = [
    ("format-version", '1.2'),
    ("date", str(datetime.datetime.now())),
    ('remark', 'namespace: MONO'),
    ('remark', 'creator: Joshua Klein <jaklein <-at-> bu.edu>'),
    ('ontology', 'MONO')
]


write_header(None, header, buff)

for key, term in sorted(cv.items(), key=lambda x: ("generic" not in x[1].definition, x[1].has_monoisotopic_mass)):
    write_term(None, term, buff)

print(buff.getvalue())
