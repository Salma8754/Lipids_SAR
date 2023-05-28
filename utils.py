
from rdkit import Chem
from rdkit.Chem import Draw
import io

# Define a function to convert a molecule string to an image

def smiles_to_img(mol_str, size=(400, 400)):
    mol = Chem.MolFromSmiles(mol_str)
    if mol is None:
        raise ValueError(f'Invalid molecule: {mol_str}')
    img = Draw.MolToImage(mol, size=size, returnPNG=False )
    return img


def smiles_to_bytes(mol_str, size=(300, 300)):
    mol = Chem.MolFromSmiles(mol_str)
    if mol is None:
        raise ValueError(f'Invalid molecule: {mol_str}')
    img = Draw.MolToImage(mol, size=size)
    image_buffer = io.BytesIO()
    img.save(image_buffer, format="JPEG")
    return img


def smarts_to_img(mol_str, size=(300, 300)):
    mol = Chem.MolFromSmarts(mol_str)
    if mol is None:
        raise ValueError(f'Invalid molecule: {mol_str}')
    img = Draw.MolToImage(mol, size=size)
    return img