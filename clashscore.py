# -*- coding: utf-8 -*-
from Bio import PDB
import numpy as np
from Bio.PDB import Superimposer
import os
import sys

def calculate_clash_score(pdb_filename):
    # Wczytaj struktur� PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_filename)

    clash_score = 0

    # Iteruj przez modele, �a�cuchy i reszty
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Sprawd� kolizje z innymi atomami w strukturze
                    for other_model in structure:
                        for other_chain in other_model:
                            for other_residue in other_chain:
                                for other_atom in other_residue:
                                    # Pomijaj por�wnywanie atom�w w tej samej reszcie
                                    if atom.id != other_atom.id:
                                        # Oblicz odleg�o�� mi�dzy atomami
                                        distance = atom - other_atom
                                        
                                        # Ustal pr�g kolizji (np. 2.0 ?ngstr�ma)
                                        clash_threshold = 2.0

                                        # Je�li odleg�o�� jest mniejsza ni� pr�g, dodaj do Clash Score
                                        if distance < clash_threshold:
                                            clash_score += 1

    return clash_score

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uzycie: python script.py nazwa_pliku.pdb")
        sys.exit(1)

    pdb_filename = sys.argv[1]
    result = calculate_clash_score(pdb_filename)
    print(f"Clash Score: {result}")

