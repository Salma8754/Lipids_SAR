import streamlit as st
from streamlit_image_select import image_select

from PIL import Image
import PIL.Image
import os

from data_fetcher import df
from fragments import link_fragments, head_fragments
from utils import smiles_to_img, smarts_to_img


base_dir = os.path.abspath(os.path.dirname(__file__))

df["num_tails"] = df["num_tails"].astype(int)


# Set the page configuration
st.set_page_config( layout='wide')

st.markdown(
    """
    <style>
    .vertical-line {
        border-left: 1px solid #ccc;
        height: 100%;
        margin: 0 10px;
    }
    </style>
    """,
    unsafe_allow_html=True
)



# Add a title
st.markdown("## <center> Structure Activity Relationship</center>", unsafe_allow_html=True)
st.write('---')



# Create a two-column layout
col1, col2, col3 = st.columns([1, 1, 1])


# Head fragment drawings
head_imgs = [smiles_to_img(hf) for hf in head_fragments]
head_img_dict = {head_fragments[i] : head_imgs[i] for i in range(len(head_imgs))}


# Tail saturation conf
img_saturated_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'SaturatedTail.png')).convert('RGB')

img_unsaturated_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'UnsaturatedTail.png')).convert('RGB')

tail_saturation_img_dict = {'saturated_tail': img_saturated_tail, 'unsaturated_tail': img_unsaturated_tail}

# Tail branches conf
img_branched_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'BranchedTail.png')).convert('RGB')

img_unbranched_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'UnbranchedTail.png')).convert('RGB')

tail_branch_img_dict = {'branched_tail': img_branched_tail, 'unbranched_tail': img_unbranched_tail}





with col1:
    st.markdown("## Head Group")
    st.write('---')
    selected_head = image_select("Head Fragment", head_imgs, use_container_width=False)
    selected_head_smiles = list(head_img_dict.keys())[list(head_img_dict.values()).index(selected_head)]

# Filtering df with head fragment
filtered_df = df[df[selected_head_smiles]]

filtered_linkers = []
for l in link_fragments:
    if sum(filtered_df[l]) > 0:
        filtered_linkers.append(l)

# Linker drawings
link_imgs = [smarts_to_img(hf) for hf in filtered_linkers]
link_img_dict = {filtered_linkers[i] : link_imgs[i] for i in range(len(link_imgs))}

with col2:
    st.markdown("## Linker")
    st.write('---')
    if len(filtered_linkers)>0:
        selected_link = image_select("Linker", link_imgs, use_container_width=False)
        selected_linker_smiles = list(link_img_dict.keys())[list(link_img_dict.values()).index(selected_link)]


# Filtering df with linker
if len(filtered_linkers)>0:
    filtered_df = filtered_df[filtered_df[selected_linker_smiles]]

# Tail saturation indicators
sat = sum(filtered_df.double_bonds_per_tail == 0)
unsat = sum(filtered_df.double_bonds_per_tail > 0)

with col3:
    st.markdown("## Tail")
    st.write('---')
    st.markdown("#### Tail Saturation")
    tail_saturation = image_select("", [img_saturated_tail, img_unsaturated_tail], captions=[f"{sat} lipids", f"{unsat} lipids"])
    tail_saturation_msg = list(tail_saturation_img_dict.keys())[list(tail_saturation_img_dict.values()).index(tail_saturation)]

# Filtering df with tail saturation
if tail_saturation_msg == 'unsaturated_tail':
    filtered_df = filtered_df[filtered_df["double_bonds_per_tail"] > 0]
else:
    filtered_df = filtered_df[filtered_df["double_bonds_per_tail"] <= 0]

# Tail branching indicators
bran = sum(filtered_df.branches_in_tails > 0)
unbran = sum(filtered_df.branches_in_tails <= 0)


with col3:
    st.markdown("#### Tail Branching")
    tail_branch = image_select('', [img_unbranched_tail, img_branched_tail], captions=[f"{unbran} lipids", f"{bran} lipids"])
    tail_branch_msg = list(tail_branch_img_dict.keys())[list(tail_branch_img_dict.values()).index(tail_branch)]


# Filtering df with tail branching
if tail_branch_msg == 'branched_tail':
    filtered_df = filtered_df[filtered_df["branches_in_tails"] > 0]
else:
    filtered_df = filtered_df[filtered_df["branches_in_tails"] <= 0]


# Get available numbers of tails
list_tail_numbers = filtered_df["num_tails"].unique()

with col3:
    num_tails = st.selectbox('#### Number of Tails', list_tail_numbers)

# Filter df with number of tails

if num_tails:
    filtered_df = filtered_df[filtered_df["num_tails"] == num_tails]


#st.dataframe(filtered_df)

# Display the dataframe in Streamlit
st.write('---')

st.markdown("## <center> Resulting lipids</center>", unsafe_allow_html=True)


# for index, row in filtered_df.iterrows():
#     st.write(f"**Lipid ID:** {row['catlipid_id']}")
#     st.write(f"**Family:** {row['CATLIPID_FAMILY_NAME']}")
#     st.write(f"**EPO 6h :** {row['max_epo_conc_mean_6hr']} ng/ml")
#     st.write(f"**EPO 24h :** {row['max_epo_conc_mean_24hr']} ng/ml")
#
#     img = Image.open(f"C:/Users/stafasca/Documents/SANOFI/App_lipids_SAR/lipid_images/lipid{index}.jpg")
#
#     # image_stream = BytesIO()
#     # row['lipid_img'].save(image_stream, format='PNG')
#     st.image(img)
#     st.write('---')



# Display the dataframe in Streamlit




st.write('---')
for index, row in filtered_df.iterrows():
    container = st.container()
    img = Image.open(os.path.join(base_dir, 'lipid_images', f'lipid{index}.jpg'))


    with container:
        st.write(f"     **Lipid ID:** {row['catlipid_id']}")
        st.write(f"     **Family:** {row['CATLIPID_FAMILY_NAME']}")
        st.write(f"     **EPO 6h :** {row['max_epo_conc_mean_6hr']} ng/ml")
        st.write(f"     **EPO 24h :** {row['max_epo_conc_mean_24hr']} ng/ml")
        st.image(img)
        st.write('---')




# Custom CSS styling for the container
container_style = """
    border: 1px solid #ccc;
    border-radius: 5px;
    padding: 10px;
    margin-bottom: 10px;
"""