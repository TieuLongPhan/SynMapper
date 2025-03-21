from rdkit import Chem


def clean_smiles(smiles: str, canonical: bool = True) -> str:
    """
    Converts SMILES to canonical form and removes atom mappings.

    Parameters:
    - smiles (str): SMILES string.
    - canonical (bool): Whether to return canonical SMILES. Defaults to True.

    Returns:
    - str: Cleaned SMILES without atom mappings.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    # Remove atom mappings
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    return Chem.MolToSmiles(mol, canonical=canonical)


def remove_atom_mapping(reaction_smiles: str, symbol: str = ">>") -> str:
    """
    Removes atom mappings from reaction SMILES.

    Parameters:
    - reaction_smiles (str): Reaction SMILES with atom mappings.
    - symbol (str): Symbol separating reactants and products (default '>>').

    Returns:
    - str: Reaction SMILES without atom mappings.
    """
    try:
        reactants, products = reaction_smiles.split(symbol)
    except ValueError:
        raise ValueError(
            "Invalid reaction SMILES format. Expected format 'reactants>>products'."
        )

    clean_reactants = clean_smiles(reactants)
    clean_products = clean_smiles(products)

    return f"{clean_reactants}{symbol}{clean_products}"
