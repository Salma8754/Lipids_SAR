import os

import pandas as pd
# from sqlalchemy import desc
# from sqlalchemy.sql import func
# import numpy as np
# import ast
#
# from data_client import session
# from schemas import EpoMerged, catlipids
#
#
#
# def get_data():
#     """
#     get lipids data from Snowflake
#     :return:
#     --------
#     a list of lipids
#     """
#     data = (
#         session.query(
#             EpoMerged.catlipid_id,
#             func.max(EpoMerged.epo_conc_mean_6hr).label('max_epo_conc_mean_6hr'),
#             func.max(EpoMerged.epo_conc_mean_6hr).label('max_epo_conc_mean_24hr'),
#             EpoMerged.canonical_smiles,
#             EpoMerged.species,
#             catlipids.catlipid_family_name,
#             catlipids.mol_wt
#          )
#         .where(EpoMerged.composition_cat == 40.0)  # primary screening
#         .where(EpoMerged.composition_peg == 1.5)  # primary screening
#         .where(EpoMerged.composition_chol == 28.5)  # primary screening
#         .where(EpoMerged.composition_help == 30)  # primary screening
#         .join(
#              catlipids,
#              catlipids.catlipid_id == EpoMerged.catlipid_id)
#         .group_by(
#             EpoMerged.catlipid_id,
#             EpoMerged.canonical_smiles,
#             EpoMerged.species,
#             catlipids.mol_wt,
#             catlipids.catlipid_family_name).order_by(desc('max_epo_conc_mean_6hr'))
#     )
#     return data.all()

df = pd.read_csv('df_fg_tail_desc.csv')






