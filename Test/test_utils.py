import unittest
from rdkit import Chem
from synmapper.utils import clean_smiles, remove_atom_mapping


class TestSmilesFunctions(unittest.TestCase):

    def test_clean_smiles_valid(self):
        self.assertEqual(clean_smiles("[CH3:1][OH:2]"), "CO")
        self.assertEqual(
            clean_smiles(Chem.CanonSmiles("C[CH](N)C(=O)O")),
            Chem.CanonSmiles("N[CH](C)C(=O)O"),
        )

    def test_clean_smiles_invalid(self):
        with self.assertRaises(ValueError):
            clean_smiles("invalid_smiles")

    def test_remove_atom_mapping_valid(self):
        reaction_smiles = "[CH3:1][OH:2]>>[CH2:1]=[O:2]"
        cleaned_reaction = remove_atom_mapping(reaction_smiles)
        self.assertEqual(cleaned_reaction, "CO>>C=O")

    def test_remove_atom_mapping_invalid_format(self):
        with self.assertRaises(ValueError):
            remove_atom_mapping("InvalidReactionSMILES")

    def test_remove_atom_mapping_custom_symbol(self):
        reaction_smiles = "[CH3:1][OH:2]>[O]>[CH2:1]=[O:2]"
        cleaned_reaction = remove_atom_mapping(reaction_smiles, symbol=">[O]>")
        self.assertEqual(cleaned_reaction, "CO>[O]>C=O")


if __name__ == "__main__":
    unittest.main()
