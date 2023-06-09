import streamlit as st
from streamlit_image_select import image_select
import numpy as np
from plotly import express as px
from PIL import Image
import PIL.Image
import os
from data_fetcher import df
from fragments import link_fragments, head_fragments
from utils import smiles_to_img, smarts_to_img, generate_image_with_text
from streamlit_plotly_events import plotly_events
from dim_reduction import df_pca, df_umap, custom_plot_pca, custom_plot_umap
from stmol import showmol
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import streamlit.components.v1 as components


base_dir = os.path.abspath(os.path.dirname(__file__))

df["num_tails"] = df["num_tails"].astype(int)


#st.write(time.time()-itime)
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



#custom_plot_umap = cp2.interactive_plot(kind='scatter')
#custom_plot_pca = cp1.interactive_plot(kind='scatter')




# Add a title
st.markdown("## <center> Structure Activity Relationship</center>", unsafe_allow_html=True)
st.write('---')



# Create a two-column layout
col1, col2, col3 = st.columns([1, 1, 1])
img_all = generate_image_with_text('All')


# Head fragment drawings
head_imgs = [img_all] + [smiles_to_img(hf) for hf in head_fragments]
head_fragments = ['ALL'] + head_fragments
head_img_dict = {head_fragments[i] : head_imgs[i] for i in range(len(head_imgs))}


# Tail saturation conf
img_saturated_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'SaturatedTail.png')).convert('RGB')

img_unsaturated_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'UnsaturatedTail.png')).convert('RGB')

tail_saturation_img_dict = {'saturated_tail': img_saturated_tail, 'unsaturated_tail': img_unsaturated_tail}

# Tail branches conf
img_branched_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'BranchedTail.png')).convert('RGB')

img_unbranched_tail = PIL.Image.open(os.path.join(base_dir, 'st_images', 'UnbranchedTail.png')).convert('RGB')

tail_branch_img_dict = {'branched_tail': img_branched_tail, 'unbranched_tail': img_unbranched_tail}



#########################################################

#                       Head Group

#########################################################


with col1:
    st.markdown("### Head Group")
    st.write('---')
    selected_head = image_select("Select Head Fragment", head_imgs, use_container_width=False)
    selected_head_smiles = list(head_img_dict.keys())[list(head_img_dict.values()).index(selected_head)]

# Filtering df with head fragment if 'ALL' is not selected
if selected_head_smiles != 'ALL':
    filtered_df = df[df[selected_head_smiles]]
else:
    filtered_df = df

filtered_linkers = []
for l in link_fragments:
    if sum(filtered_df[l]) > 0:
        filtered_linkers.append(l)


# Linker drawings



link_imgs = [img_all] + [smarts_to_img(hf) for hf in filtered_linkers]
filtered_linkers = ['ALL'] + filtered_linkers
link_img_dict = {filtered_linkers[i] : link_imgs[i] for i in range(len(link_imgs))}


#########################################################

#                       Linkers

#########################################################

with col2:
    st.markdown("### Linker")
    st.write('---')
    if len(filtered_linkers)>0:
        selected_link = image_select("Select Linker Pattern", link_imgs, use_container_width=False)
        selected_linker_smiles = list(link_img_dict.keys())[list(link_img_dict.values()).index(selected_link)]


# Filtering df with linker
if selected_linker_smiles != 'ALL':
    print('im here')
    filtered_df = filtered_df[filtered_df[selected_linker_smiles]]

possible_deg_sat = np.sort(filtered_df.double_bonds_per_tail.unique())

#########################################################

#                       Tails

#########################################################

with col3:
    st.markdown("### Tail")
    st.write('---')
    st.markdown("#### Tail Saturation")

    selected_number_double_bonds = st.selectbox('Number of double bonds per tail (total number of double bonds normalized by number of tails)', possible_deg_sat)

if selected_number_double_bonds is not None:
    filtered_df = filtered_df[filtered_df["double_bonds_per_tail"] == selected_number_double_bonds]

# Tail saturation indicators --> not very useful
#sat = sum(filtered_df.double_bonds_per_tail == 0)
#unsat = sum(filtered_df.double_bonds_per_tail > 0)

# with col3: # not very useful
#     tail_saturation = image_select("", [img_saturated_tail, img_unsaturated_tail], captions=[f"{sat} lipids", f"{unsat} lipids"])
#     tail_saturation_msg = list(tail_saturation_img_dict.keys())[list(tail_saturation_img_dict.values()).index(tail_saturation)]
#
# # Filtering df with tail saturation
# if tail_saturation_msg == 'unsaturated_tail':
#     filtered_df = filtered_df[filtered_df["double_bonds_per_tail"] > 0]
# else:
#     filtered_df = filtered_df[filtered_df["double_bonds_per_tail"] <= 0]

# Tail branching indicators
bran = sum(filtered_df.branches_in_tails > 0)
unbran = sum(filtered_df.branches_in_tails <= 0)


#########################################################

#                       Tails

#########################################################

with col3:
    st.markdown("#### Tail Branching")
    tail_branch = image_select('', [img_unbranched_tail, img_branched_tail], captions=[f"{unbran} lipids", f"{bran} lipids"], use_container_width=False)
    tail_branch_msg = list(tail_branch_img_dict.keys())[list(tail_branch_img_dict.values()).index(tail_branch)]


# Filtering df with tail branching
if tail_branch_msg == 'branched_tail':
    filtered_df = filtered_df[filtered_df["branches_in_tails"] > 0]
else:
    filtered_df = filtered_df[filtered_df["branches_in_tails"] <= 0]


# Get available numbers of tails
list_tail_numbers = np.sort(filtered_df["num_tails"].unique())

with col3:
    st.markdown("#### Number of tails")
    num_tails = st.selectbox('', list_tail_numbers)

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





#scatter = px.scatter(filtered_df, x=filtered_df.max_epo_conc_mean_6hr, y=filtered_df.max_epo_conc_mean_6hr, color=filtered_df.max_epo_conc_mean_6hr, color_continuous_scale='viridis', opacity=0.7, title="Lipids Space")
scatter1 = px.scatter(df_pca, x=df_pca.columns[0], y=df_pca.columns[1], color=df_pca["target"], color_continuous_scale='viridis', opacity=0.7, title="Lipids Space")
#scatter2 = px.scatter(df_tsne, x=df_tsne.columns[0], y=df_tsne.columns[1], color=df_tsne["target"], color_continuous_scale='viridis', opacity=0.7, title="Lipids Space")
scatter3 = px.scatter(df_umap, x=df_umap.columns[0], y=df_umap.columns[1],
                      color=df_umap["target"], color_continuous_scale='viridis',
                      opacity=0.7, title="Lipids Space"                      )

#scatter3.update_traces(customdata=lip_imgs)

scatter1.update_layout(height=400, width=500)
#scatter2.update_layout(height=400, width=500)
scatter3.update_layout(height=400, width=500)














#########################################################

#                       SIDEBAR

#########################################################

st.sidebar.write('**Visualization Parameters**')

sim_type = st.sidebar.radio(
    "Which similarity type do you want to use?",
    ('tailored', 'structural'),
    help='Use tailored when you have a target value. Use structural to plot your molecules based on structure only.')

dim_red_algo = st.sidebar.radio(
    "Which algorithm you want to use?",
    ('PCA', 'UMAP'),
    help='UMAP is non linear resulting in better clusters. PCA is linear resulting in a more global view.')

plot_type = st.sidebar.radio(
    "Which plot type do you want to display?",
    ('scatter', 'hex'),
    help='Visualize a scatter plot or an hexagonal plot.')

#########################################################

#                       Visualization

#########################################################


col4, col5 = st.columns((1, 3))

hvar = """ 
        <script> 
                var elements = window.parent.document.querySelectorAll('.streamlit-expanderHeader');
                elements[0].style.color = 'rgba(83, 36, 118, 1)';
                elements[0].style.fontSize = 'x-large';
        </script>
    """

data_scatter_plot = st.expander("Visualize the Lipids Space", expanded=False)

with data_scatter_plot:
    if dim_red_algo == 'PCA':

        st.bokeh_chart(custom_plot_pca, use_container_width=True)
        #selected_points1 = plotly_events(scatter1, select_event=True, click_event=True, hover_event=False)
    else:
        st.bokeh_chart(custom_plot_umap, use_container_width=True)
        #selected_points2 = plotly_events(scatter2, select_event=True, click_event=True, hover_event=False)
        #selected_points3 = plotly_events(scatter3, select_event=True, click_event=True, hover_event=False)



data_lipids_plot = st.expander("Visualize individual Lipids", expanded=False)

with col4:
    with data_lipids_plot:
        st.write('---')
        for index, row in filtered_df.iterrows():
            container = st.container()
            img = Image.open('lipid_images/' + f'lipid{index}.jpg')

            with container:
                st.write(f"     **Lipid ID:** {row['catlipid_id']}")
                st.write(f"     **Family:** {row['CATLIPID_FAMILY_NAME']}")
                st.write(f"     **EPO 6h :** {row['max_epo_conc_mean_6hr']} ng/ml")
                st.write(f"     **EPO 24h :** {row['max_epo_conc_mean_24hr']} ng/ml")
                st.image(img)
                st.write('---')



hvar = """ 
    <script> 
            var elements = window.parent.document.querySelectorAll('.streamlit-expanderHeader');
            elements[0].style.color = 'rgba(83, 36, 118, 1)';
            elements[0].style.fontSize = 'x-large';
            elements[1].style.color = 'rgba(83, 36, 118, 1)';
            elements[1].style.fontSize = 'x-large';
    </script>
"""

components.html(hvar, height=0, width=0)

########### 3D mol plot ##############

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    #mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    xyzview = py3Dmol.view()#(width=400,height=400)
    xyzview.addModel(xyz,'mol')
    xyzview.setStyle({'stick': {}})

    # Add surface representation with contour --> crazy
    #xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.8, 'color': 'spectrum', 'contour': True})

    xyzview.setBackgroundColor('white')
    #xyzview.zoomTo()

    showmol(xyzview, height=500, width=500)


#compound_smiles=st.text_input('SMILES please','CC')
#blk=makeblock(compound_smiles)
#render_mol(blk)

########### end 3D mol plot ##############

