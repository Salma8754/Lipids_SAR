import os

import streamlit as st
from streamlit_image_select import image_select
import streamlit.components.v1 as components

import numpy as np
from plotly import express as px
from PIL import Image
import PIL.Image
from PIL import Image, ImageDraw, ImageFont

from fragments import link_fragments, head_fragments

from utils import smiles_to_img, smarts_to_img, makeblock, render_mol, generate_image_with_text
from dim_reduction import df_pca, df_umap, custom_plot_pca, custom_plot_umap
from data_fetcher import df




base_dir = os.path.abspath(os.path.dirname(__file__))

df["num_tails"] = df["num_tails"].astype(int)



#########################################################

#                       Layout

#########################################################


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
st.markdown("## <center> Exploring Lipids Space - Structure Activity Relationship</center>", unsafe_allow_html=True)
st.markdown("##### <center> Filter Lipids by Head-Linker-Tail properties</center>", unsafe_allow_html=True)

st.markdown("## ", unsafe_allow_html=True)

st.write('---')

st.markdown("## ", unsafe_allow_html=True)

#st.write('---')


# Create a two-column layout
col1, col2, col3 = st.columns([1, 1, 1])




img_all = generate_image_with_text('Unspecified')


#########################################################

#                       SIDEBAR

#########################################################

st.sidebar.write('**Visualization Parameters**')

dim_red_algo = st.sidebar.radio(
            "Algorithm used for the scatterplot",
            ('PCA', 'UMAP'),
            help='UMAP is non linear resulting in better clusters. PCA is linear resulting in a more global view.')

dim_visual_mol = st.sidebar.radio(
            "Lipid visualization type",
            ('2D', '3D', '3D surface representation'))





#########################################################

#               Head-Linker-Tail config

#########################################################


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





#########################################################

#                       Linkers

#########################################################

# Linker drawings
link_imgs = [img_all] + [smarts_to_img(hf) for hf in filtered_linkers]
filtered_linkers = ['ALL'] + filtered_linkers
link_img_dict = {filtered_linkers[i] : link_imgs[i] for i in range(len(link_imgs))}

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


#########################################################

#                       Tails

#########################################################

# Tail sat config
possible_deg_sat = np.sort(filtered_df.double_bonds_per_tail.unique())


with col3:
    st.markdown("### Tail")
    st.write('---')
    st.markdown("#### Tail Saturation")

    selected_number_double_bonds = st.selectbox('Number of double bonds per tail (total number of double bonds normalized by number of tails)', ['Unspecified']+list(possible_deg_sat))

if selected_number_double_bonds is not None and selected_number_double_bonds != 'Unspecified':
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


# Tail branches config
possible_num_branch = np.sort(filtered_df.branches_in_tails.unique())

# Tail branching indicators
bran = sum(filtered_df.branches_in_tails > 0)
unbran = sum(filtered_df.branches_in_tails <= 0)


with col3:
    st.markdown("#### Branched-Tail(s)")
    selected_number_branches = st.selectbox('Whether there are branches in the lipid tails', ['Unspecified', 'Branched-Tails', 'Unbranched-Tails'])

    #tail_branch = image_select('', [img_unbranched_tail, img_branched_tail], captions=[f"{unbran} lipids", f"{bran} lipids"], use_container_width=False)
    #tail_branch_msg = list(tail_branch_img_dict.keys())[list(tail_branch_img_dict.values()).index(tail_branch)]


# Filtering df with tail branching
#if tail_branch_msg == 'branched_tail':
    #filtered_df = filtered_df[filtered_df["branches_in_tails"] > 0]
#else:
    #filtered_df = filtered_df[filtered_df["branches_in_tails"] <= 0]

if selected_number_branches is not None and selected_number_branches == 'Branched-Tails':
    filtered_df = filtered_df[filtered_df["branches_in_tails"] == 1]
elif selected_number_branches is not None and selected_number_branches == 'Unbranched-Tails':
    filtered_df = filtered_df[filtered_df["branches_in_tails"] == 0]



# Get available numbers of tails
list_tail_numbers = np.sort(filtered_df["num_tails"].unique())

with col3:
    st.markdown("#### Number of tails")
    num_tails = st.selectbox('', ['Unspecified']+list(list_tail_numbers))

# Filter df with number of tails
if num_tails and num_tails != 'Unspecified':
    filtered_df = filtered_df[filtered_df["num_tails"] == num_tails]




# Filter dataframe with descending epo

#filtered_df = filtered_df.sort_values('max_epo_conc_mean_6hr')


#########################################################

#                   Results layout

#########################################################


# Display the dataframe in Streamlit
st.write('---')
#st.markdown(f"## <center> Resulting lipids ({len(filtered_df)}) </center>", unsafe_allow_html=True)



#########################################################

#                       Display Individual Lipids

#########################################################


data_lipids_plot = st.expander(f"Resulting Lipids ({len(filtered_df)})", expanded=False)
with data_lipids_plot:
    st.write('---')
    for index, row in filtered_df.head(10).iterrows():
        container = st.container()
        img = Image.open('lipid_images/' + f'lipid{index}.jpg')

        with container:
            st.write(f"     **Lipid ID:** {row['catlipid_id']}")
            st.write(f"     **Family:** {row['CATLIPID_FAMILY_NAME']}")
            st.write(f"     **EPO 6h :** {row['max_epo_conc_mean_6hr']} ng/ml")
            st.write(f"     **EPO 24h :** {row['max_epo_conc_mean_24hr']} ng/ml")
            if dim_visual_mol == '2D':
                st.image(img)
            elif dim_visual_mol == '3D':
                compound_smiles = row['canonical_smiles']
                blk = makeblock(compound_smiles)
                render_mol(blk, surf=False)
            else:
                compound_smiles = row['canonical_smiles']
                blk = makeblock(compound_smiles)
                render_mol(blk, surf=True)
            st.write('---')






#########################################################

#                   Display Scatter  Plots

#########################################################

filtered_df_pca = df_pca[df_pca.index.isin(filtered_df.index)]
filtered_df_umap = df_umap[df_umap.index.isin(filtered_df.index)]

#scatter = px.scatter(filtered_df, x=filtered_df.max_epo_conc_mean_6hr, y=filtered_df.max_epo_conc_mean_6hr, color=filtered_df.max_epo_conc_mean_6hr, color_continuous_scale='viridis', opacity=0.7, title="Lipids Space")
scatter1 = px.scatter(filtered_df_pca, x=filtered_df_pca.columns[0], y=filtered_df_pca.columns[1], color=filtered_df_pca["target"], color_continuous_scale='viridis', opacity=0.7, title="Lipids Space")
#scatter2 = px.scatter(df_tsne, x=df_tsne.columns[0], y=df_tsne.columns[1], color=df_tsne["target"], color_continuous_scale='viridis', opacity=0.7, title="Lipids Space")
scatter3 = px.scatter(filtered_df_umap, x=filtered_df_umap.columns[0], y=filtered_df_umap.columns[1],
                      color=filtered_df_umap["target"], color_continuous_scale='viridis',
                      opacity=0.7, title="Lipids Space"                      )


scatter1.update_layout(height=400, width=500)
#scatter2.update_layout(height=400, width=500)
scatter3.update_layout(height=400, width=500)



data_scatter_plot = st.expander("Visualize the Lipids Space", expanded=False)

with data_scatter_plot:
    st.markdown("<h6 style='text-align: center;'>EPO concentration 6hr (ng/ml)</h6>", unsafe_allow_html=True)
    if dim_red_algo == 'PCA':
        st.bokeh_chart(custom_plot_pca, use_container_width=True)
        #selected_points1 = plotly_events(scatter1, select_event=True, click_event=True, hover_event=False)
    else:
        st.bokeh_chart(custom_plot_umap, use_container_width=True)
        #selected_points2 = plotly_events(scatter2, select_event=True, click_event=True, hover_event=False)
        #selected_points3 = plotly_events(scatter3, select_event=True, click_event=True, hover_event=False)




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


#########################################################

#                    Footnotes

#########################################################


st.write('---')
st.write('##### References')
st.write('Interactive app:       https://akshay.bio/blog/interactive-browser')
st.write('Umap and PCA visuals:   https://github.com/mcsorkun/ChemPlot-web')

st.write('**App in progress, more notes and features to come.**')





