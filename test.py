import pydoStuff
from pydoStuff import *

# sanity check for multiple copies of libMØD
modValue = mod.magicLibraryValue()
ourValue = pydoStuff.magicLibraryValue()
if modValue != ourValue:
	print("mod =", modValue)
	print("our =", ourValue)
	raise Exception("Magic values differ! I.e., more than one instance of libmod has been loaded.")
# end of sanity check

# make doStuff a bit more friendly to use
_doStuff_orig = pydoStuff.doStuff
def _doStuff(educts, products, doChemistryCheck=True, k=1, c=6):
	return mod._unwrap(_doStuff_orig(
		mod._wrap(mod.libpymod._VecGraph, educts),
		mod._wrap(mod.libpymod._VecGraph, products),
		doChemistryCheck,
		k,
		c))

pydoStuff.doStuff = _doStuff
doStuff = _doStuff
# end of friendlyfier code

# educts = [smiles("OCC=O")]
# products = [smiles("OC=CO")]

educts = [smiles("C=C"), smiles("C=C")]
products = [smiles("C1CCC1")]

# educts = [smiles("O"), smiles("Cl"), smiles("CC(=O)OCC")]
# products = [smiles("Cl"), smiles("OCC"), smiles("CC(=O)O")]

educts = [smiles("C1C(O)CC(O)C(O)C1")]
products = [smiles("C=CO"), smiles("C=CO"), smiles("C=CO")]

res = doStuff(educts, products, k=2, c=6)
for a in res:
	a.print()
