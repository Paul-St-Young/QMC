import numpy as np
import re

class Molecule:
    
    def __init__(self,name="molecule"):
        self.name = name
        self.formula = None
        self.atom_list = None
        self.known_atoms = {"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6\
            ,"N":7,"O":8,"F":9,"Ne":10,"Cl":17}
    # end def __init__
    
    def read_formula(self,formula):
        self.formula = formula
    # end def read_formula  
    
    def parse_formula(self,formula):
        
        # Parsing chemical formula by NPE on StackOverflow 
        first_pass = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        
        atom_list = []
        for item in first_pass:
            name,number = item
            if number == '':
                number = 1
            else:
                number = int(number)
            # end if
            atom_list.append((name,number))
        # end for item
        
        if self.formula == None:
            self.formula = formula
        # end if
        self.atom_list = atom_list
        
        return atom_list
                
    # end def parse_formula
    
    def pretty_formula(self,formula):
        atom_list = self.parse_formula(formula)
        pretty = ""
        for item in atom_list:
            name,num = item
            if num != 1:
                pretty += name + "$_" + str(num) + "$"
            else:
                pretty += name
            # end if
        # end for item
        return pretty
    # end def pretty_formula
        
    
    def num_electrons(self,formula=None):
        if formula == None:
            formula = self.formula
            atom_list = self.atom_list
        else:
            atom_list = self.parse_formula(formula)
        # end if
        
        nel = 0
        for item in atom_list:
            name,num = item
            nel += self.known_atoms[name] * num
        # end for item
        
        return nel
    # end def num_electrons
    
# end class
