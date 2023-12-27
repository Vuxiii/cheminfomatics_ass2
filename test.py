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
# end of friendlyfier code

# educts = [smiles("OCC=O")]
# products = [smiles("OC=CO")]

# educts = [smiles("O"), smiles("Cl"), smiles("CC(=O)OCC")]
# products = [smiles("Cl"), smiles("OCC"), smiles("CC(=O)O")]

# educts = [smiles("C1C(O)CC(O)C(O)C1")]
# products = [smiles("C=CO"), smiles("C=CO"), smiles("C=CO")]

# res = doStuff(educts, products, k=2, c=6)
# for a in res:
# 	a.print()

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
	
# Wrap below code inside an option which the user can run when being prompted for which test to run

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

options[18].func()
# Prompt the user for which of the above options to run. Repeat this untill all the options have been run. from std out grep for "t: " and store the number in t and grep for "n: " and store the number in n.



# doneOptions = []

# while len(options) > 0:
#     print("Choose one of the following options:")
#     for i in range(len(options)):
#         print(i, ":", options[i].educts, "->", options[i].products, "k =", options[i].k, "c =", options[i].c)
        
#     print("-"*80)

# # Print the already done options
#     for i in range(len(doneOptions)):
#         print(i, ":", doneOptions[i].educts, "->", doneOptions[i].products, "k =", doneOptions[i].k, "c =", doneOptions[i].c, "t =", doneOptions[i].t, "n =", doneOptions[i].n, "s =", doneOptions[i].s)

#     print("Or type 'done' to exit.")
#     choice = input()
#     if choice == "done":
#         break
#     try:
#         choice = int(choice)
#     except:
#         print("Invalid choice")
#         continue
#     if choice >= len(options):
#         print("Invalid choice")
#         continue
#     option = options[choice]
#     print("Running option", choice, ":", option.educts, "->", option.products, "k =", option.k, "c =", option.c)
        
#     original_out = sys.stdout
#     sys.stdout = StringIO()
        
#     option.func()

#     lines = sys.stdout.getvalue().splitlines()

#     sys.stdout = original_out

#     for line in lines:
#         if line.startswith("t: "):
#             option.t = int(line[3:])
#         elif line.startswith("n: "):
#             option.n = int(line[3:])
#         elif line.startswith("sys"):
#             option.s = str(line[4:])


#     doneOptions.append(option)
#     options.remove(option)