from io import StringIO
import pydoStuff
from pydoStuff import *
import sys

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


class Option:
    def __init__(self, educts:list, products:list, k:int, c:int):
        self.educts = [smiles(educt) for educt in educts]
        self.products = [smiles(product) for product in products]
        self.k = k
        self.c = c
        
    def func(self): 
          res = doStuff(self.educts, self.products, k=self.k, c=self.c)
          for a in res:
                a.print()

options = [ 
    Option(["C=C", "C=C"], ["C1CCC1"], 0, 4),
    Option(["C=C", "C=C"], ["C1CCC1"], 1, 4),
    Option(["O", "Cl", "CC(=O)OCC"], ["Cl", "OCC", "CC(=O)O"], 0, 6),
    Option(["O", "Cl", "CC(=O)OCC"], ["Cl", "OCC", "CC(=O)O"], 1, 6),
    Option(["O", "Cl", "CC(=O)OCC"], ["Cl", "OCC", "CC(=O)O"], 2, 6),
    Option(["C1C(O)CC(O)C(O)C1"], ["C=CO", "C=CO", "C=CO"], 0, 6),
    Option(["C1C(O)CC(O)C(O)C1"], ["C=CO", "C=CO", "C=CO"], 1, 6),
    Option(["CC=CC=CC", "OC1C=CC=CC=1"], ["C=CC=CC=C", "OC(=C)C=CC=C"], 0, 6),
    Option(["CC=CC=CC", "OC1C=CC=CC=1"], ["C=CC=CC=C", "OC(=C)C=CC=C"], 1, 6),
    Option(["CC=CC=CC", "OC1C=CC=CC=1"], ["C=CC=CC=C", "OC(=C)C=CC=C"], 2, 6),
    Option(["CC", "OC1C=CC=CC=1"], ["C=C", "OC(=C)C=CC=C"], 0, 6),
    Option(["CC", "OC1C=CC=CC=1"], ["C=C", "OC(=C)C=CC=C"], 1, 6),
    Option(["CC", "OC1C=CC=CC=1"], ["C=C", "OC(=C)C=CC=C"], 2, 6),
    Option(["OP(=O)(O)OP(=O)(O)O", "O"], ["O=P(O)(O)O", "O=P(O)(O)O"], 0, 4),
    Option(["OP(=O)(O)OP(=O)(O)O", "O"], ["O=P(O)(O)O", "O=P(O)(O)O"], 1, 4),
    Option(["OP(=O)(O)OP(=O)(O)O", "O"], ["O=P(O)(O)O", "O=P(O)(O)O"], 0, 6),
    Option(["OP(=O)(O)OP(=O)(O)O", "O"], ["O=P(O)(O)O", "O=P(O)(O)O"], 1, 6),
    Option(["C#N", "C#N"], ["N=CC#N"], 0, 4),
    Option(["C#N", "C#N"], ["N=CC#N"], 1, 4)
]

def example(k):
    educt = [ smiles(g) for g in ["O","CC(=O)OCC"] ]
    product = [ smiles(g) for g in ["OCC", "CC(=O)O"] ]
    inputRules = doStuff(educt, product, k=k, c=6)
    

    dg = DG(graphDatabase=educt)
    dg.build().execute(
        addSubset(educt)
        >> repeat[1](inputRules)
    )
    dg.print()
    # for a in dg.products: a.print()
