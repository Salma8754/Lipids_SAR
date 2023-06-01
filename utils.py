
from rdkit import Chem
from rdkit.Chem import Draw
import io
from PIL import Image, ImageDraw, ImageFont

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


def generate_image_with_text(text, font_size=24, text_color=(0, 0, 0), background_color=(255, 255, 255),
                             image_size=(100, 100), font_path=None):
    # Create a blank image
    image = Image.new('RGB', image_size, background_color)

    # Create a drawing object
    draw = ImageDraw.Draw(image)

    # Set the font
    if font_path:
        font = ImageFont.truetype(font_path, font_size)
    else:
        font = ImageFont.load_default()

    # Calculate the position to center the text
    text_width, text_height = draw.textsize(text, font=font)
    x = (image_size[0] - text_width) // 2
    y = (image_size[1] - text_height) // 2

    # Draw the text on the image
    draw.text((x, y), text, font=font, fill=text_color)

    return image




